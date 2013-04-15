// File: cuint.cu
// Author:Tom Ostler
// Last-modified: 15 Apr 2013 12:43:40
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
    //uniform temperature with interaction matrix
    __global__ void CHeun1(int N,double T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T)*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			const double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
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
    //on-site temperature with interaction matrix
    __global__ void CHeun1(int N,double *T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T[i])*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			const double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
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
    //using neighbour list to calculate exchange field, uniform temperature
    __global__ void CHeun1(int N,double T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn,int *Cxadj,int *Cadjncy,float *CJDiag)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T)*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            for(unsigned int n = Cxadj[i] ; n < Cxadj[i+1] ; n++)
            {
                unsigned int neigh=Cadjncy[n];
                h[0]+=CJDiag[3*neigh]*  Cspin[3*neigh];
                h[1]+=CJDiag[3*neigh+1]*Cspin[3*neigh+1];
                h[2]+=CJDiag[3*neigh+2]*Cspin[3*neigh+2];
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
    //using neighbour list to calculate exchange field, on-site temperature
    __global__ void CHeun1(int N,double *T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn,int *Cxadj,int *Cadjncy,float *CJDiag)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T[i])*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            for(unsigned int n = Cxadj[i] ; n < Cxadj[i+1] ; n++)
            {
                unsigned int neigh=Cadjncy[n];
                h[0]+=CJDiag[3*neigh]*  Cspin[3*neigh];
                h[1]+=CJDiag[3*neigh+1]*Cspin[3*neigh+1];
                h[2]+=CJDiag[3*neigh+2]*Cspin[3*neigh+2];
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
    //uniform temperature with interaction matrix
    __global__ void CHeun2(int N,double T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T)*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			const double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cespin[3*i],Cespin[3*i+1],Cespin[3*i+2]};
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
    //on-site temperature with interaction matrix
    __global__ void CHeun2(int N,double *T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T[i])*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			const double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
			const double s[3]={Cespin[3*i],Cespin[3*i+1],Cespin[3*i+2]};
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
    //uniform temp & CSR neighbourlist
    __global__ void CHeun2(int N,double T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn,int *Cxadj,int *Cadjncy,float *CJDiag)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T)*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
            for(unsigned int n = Cxadj[i] ; n < Cxadj[i+1] ; n++)
            {
                unsigned int neigh=Cadjncy[n];
                h[0]+=CJDiag[3*neigh]*  Cespin[3*neigh];
                h[1]+=CJDiag[3*neigh+1]*Cespin[3*neigh+1];
                h[2]+=CJDiag[3*neigh+2]*Cespin[3*neigh+2];
            }
			const double s[3]={Cespin[3*i],Cespin[3*i+1],Cespin[3*i+2]};
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
    //on-site temp & CSR neighbourlist
    __global__ void CHeun2(int N,double *T,double sigma,double llgpf,double lambda,double rdt,double appliedx,double appliedy,double appliedz,float *CH,double *Cspin,double *Cespin,float *Crand,double *Cfn,int *Cxadj,int *Cadjncy,float *CJDiag)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			//The prefactor for the thermal term
			const double TP=sqrt(T[i])*sigma;
			const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
			double h[3]={double(CH[3*i])+lrn[0]+appliedx,double(CH[3*i+1])+lrn[1]+appliedy,double(CH[3*i+2])+lrn[2]+appliedz};
            for(unsigned int n = Cxadj[i] ; n < Cxadj[i+1] ; n++)
            {
                unsigned int neigh=Cadjncy[n];
                h[0]+=CJDiag[3*neigh]*  Cespin[3*neigh];
                h[1]+=CJDiag[3*neigh+1]*Cespin[3*neigh+1];
                h[2]+=CJDiag[3*neigh+2]*Cespin[3*neigh+2];
            }
			const double s[3]={Cespin[3*i],Cespin[3*i+1],Cespin[3*i+2]};
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
