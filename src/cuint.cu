// File: cuint.cu
// Author:Tom Ostler
// Last-modified: 23 Aug 2022 15:19:32
#include "../inc/cufields.h"
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/util.h"
#include "../inc/cuint.h"
#include "../inc/fields.h"
#include "../inc/defines.h"
#include "../inc/cudadefs.h"
#include "../inc/anis.h"
#include "../inc/llg.h"
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
namespace cuint
{
    //Reduced timestep
    __constant__ double Crdt;

    void copyConstData()
    {
        FIXOUT(config::Info,"Copying const data with cuint scope to card:" << std::flush);
        cudaMemcpyToSymbol(*(&Crdt),&llg::rdt,sizeof(double));
        config::Info << "Done" << std::endl;
    }
    __global__ void CHeun1(int N,double T,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //The prefactor for the thermal term
            const double TP=sqrt(T)*Csigma[i];
            const double llgpf=Cllgpf[i];
            const double lambda=Clambda[i];
            const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};

            double h[3]={double(CH[3*i])+CHstg[3*i]+lrn[0]+appliedx,double(CH[3*i+1])+CHstg[3*i+1]+lrn[1]+appliedy,double(CH[3*i+2])+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            //calculate the field arising from the first order uniaxial anisotropy
            const double k1udir[3]={Ck1udir[3*i],Ck1udir[3*i+1],Ck1udir[3*i+2]};
            const double k1u=Ck1u[i];
            const double sdn = s[0]*k1udir[0] + s[1]*k1udir[1] + s[2]*k1udir[2];

            h[0]+=(k1u*sdn*k1udir[0]);
            h[1]+=(k1u*sdn*k1udir[1]);
            h[2]+=(k1u*sdn*k1udir[2]);
            //printf("%4.10f\t%4.10f\t%4.10f\n",h[0],h[1],h[2]);

            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            double lfn[3]={0,0,0};
            double es[3]={0,0,0};
            double mods=0.0;
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lfn[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
                Cfn[3*i+j]=lfn[j];
                es[j]=s[j]+lfn[j]*Crdt;
                mods+=es[j]*es[j];
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

    __global__ void CHeun2(int N,double T,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CDetField,double *CHstg)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //The prefactor for the thermal term
            const double TP=sqrt(T)*Csigma[i];
            const double llgpf=Cllgpf[i];
            const double lambda=Clambda[i];
            const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
            double h[3]={double(CH[3*i])+CHstg[3*i]+lrn[0]+appliedx,double(CH[3*i+1])+CHstg[3*i+1]+lrn[1]+appliedy,double(CH[3*i+2])+CHstg[3*i+2]+lrn[2]+appliedz};
/*            if(i==1 || i==0)
            {
                printf("spin=%d\t%4.6f\t%4.6f\t%4.6f\n",i,h[0],h[1],h[2]);
            }*/

            const double s[3]={Cespin[3*i],Cespin[3*i+1],Cespin[3*i+2]};
            //calculate the field arising from the first order uniaxial anisotropy
            const double k1udir[3]={Ck1udir[3*i],Ck1udir[3*i+1],Ck1udir[3*i+2]};
            const double k1u=Ck1u[i];
            const double sdn = s[0]*k1udir[0] + s[1]*k1udir[1] + s[2]*k1udir[2];

            h[0]+=(k1u*sdn*k1udir[0]);
            h[1]+=(k1u*sdn*k1udir[1]);
            h[2]+=(k1u*sdn*k1udir[2]);
            double ps[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};
            const double fn[3]={Cfn[3*i],Cfn[3*i+1],Cfn[3*i+2]};
            double fnp1[3]={0,0,0};
            double mods=0.0;
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                fnp1[j]=llgpf*(sxh[j]+lambda*sxsxh[j]);
                ps[j]+=(0.5*(fn[j]+fnp1[j])*Crdt);
                mods+=ps[j]*ps[j];
            }
            const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                Cspin[3*i+j]=ps[j]*nf;
                CDetField[3*i+j]=h[j]-lrn[j];
                CH[3*i+j]=0.0;
            }
            //zero the CH array here


        }
    }
}
