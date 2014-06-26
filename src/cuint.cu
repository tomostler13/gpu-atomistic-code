// File: cuint.cu
// Author:Tom Ostler
// Last-modified: 26 Jun 2014 11:45:09
#include "../inc/cufields.h"
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/util.h"
#include "../inc/cuint.h"
#include "../inc/fields.h"
#include "../inc/defines.h"
#include "../inc/cudadefs.h"
#include "../inc/anis.h"
//Cuda headers
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cufft.h>
//Library headers
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#define MAXNFOU 5
namespace cuint
{
    __constant__ int Cnfou=0;
    __constant__ double CUK[MAXNFOU] = {0,0,0,0,0};
    //Will need to do a 2D array lookup on this
    __constant__ double CUKD[MAXNFOU*3] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    void copyConstData()
    {
        FIXOUT(config::Info,"Copying const data to card:" << std::flush);
        cudaMemcpyToSymbol(*(&Cnfou),&anis::nfou,sizeof(int));
        cudaMemcpyToSymbol(CUK,anis::FirstOrderUniaxK.ptr(),anis::FirstOrderUniaxK.size()*sizeof(double));
        cudaMemcpyToSymbol(CUKD,anis::FirstOrderUniaxDir.ptr(),anis::FirstOrderUniaxDir.size()*sizeof(double));
        config::Info << "Done" << std::endl;
    }
    __global__ void CHeun1(int N,double T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            /*for(unsigned int k=0;k<Cnfou;k++)
            {
                printf("Uniaxial anisotropy constant %d of %d:\t%4.10f\tin direction\t(%4.2f,%4.2f,%4.2f)\n",k,Cnfou,CUK[k],CUKD[k*3],CUKD[k*3+1],CUKD[k*3+2]);
            }*/
            //The prefactor for the thermal term
			const double TP=sqrt(T)*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};

            printf("%f\t%f\t%f\n",s[0],s[1],s[2]);
            for(unsigned int k = 0 ; k < Cnfou ; k++)
            {
                double sdotn=s[0]*CUKD[k*3]*CUK[k]+s[1]*CUKD[k*3+1]*CUK[k]+s[2]*CUKD[k*3+2]*CUK[k];
                h[0]+=sdotn*CUKD[k*3];
                h[1]+=sdotn*CUKD[k*3+1];
                h[2]+=sdotn*CUKD[k*3+2];
            }

			const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
			const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

			double lfn[3]={0,0,0};
			double es[3]={0,0,0};
			double mods=0.0;
			for(unsigned int j = 0 ; j < 3 ; j++)
			{
				lfn[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
				Cfn[3*i+j]=lfn[j];
				es[j]=s[j]+lfn[j]*rdt;
				mods+=s[j]*s[j];
			}
			//calculate one over the square root of the spin modulus
			const double nf=rsqrt(mods);
			for(unsigned int j = 0 ; j < 3 ; j++)
			{
				//set the euler spin value and normalize
				Cespin[3*i+j]=es[j]*nf;
			}
        }
    }

    __global__ void CHeun2(int N,double T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T)*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cespin[3*i],Cespin[3*i+1],Cespin[3*i+2]};
            for(unsigned int k = 0 ; k < Cnfou ; k++)
            {
                double sdotn=s[0]*CUKD[k*3]*CUK[k]+s[1]*CUKD[k*3+1]*CUK[k]+s[2]*CUKD[k*3+2]*CUK[k];
                h[0]+=sdotn*CUKD[k*3];
                h[1]+=sdotn*CUKD[k*3+1];
                h[2]+=sdotn*CUKD[k*3+2];
            }
			double ps[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
			const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
			const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};
			const double fn[3]={Cfn[3*i],Cfn[3*i+1],Cfn[3*i+2]};
            double fnp1[3]={0,0,0};
			double mods=0.0;
			for(unsigned int j = 0 ; j < 3 ; j++)
			{
				fnp1[j]=llgpf*(sxh[j]+lambda*sxsxh[j]);
				ps[j]+=(0.5*(fn[j]+fnp1[j])*rdt);
				mods+=ps[j]*ps[j];
			}
			const double nf=rsqrt(mods);
			for(unsigned int j = 0 ; j < 3 ; j++)
			{
				Cspin[3*i+j]=ps[j]*nf;
			}

        }
    }
}
