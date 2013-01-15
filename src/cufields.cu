// File: cufields.cu
// Author:Tom Ostler
// Last-modified: 03 Jan 2013 13:08:59
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
//cuda headers
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cufft.h>
//local cuda headers
#include "../inc/cufields.h"

#ifdef DEBUG
#define CUDA_CALL(x) do { if((x) != cudaSuccess) {\
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CUDA_CALL(x) (x)
#endif /*DEBUG*/
namespace cufields
{
    //perform the convolution in Fourier space
    __global__ void CFConv(int N,
                           cufftComplex *CCNxx,
                           cufftComplex *CCNxy,
                           cufftComplex *CCNxz,
                           cufftComplex *CCNyx,
                           cufftComplex *CCNyy,
                           cufftComplex *CCNyz,
                           cufftComplex *CCNzx,
                           cufftComplex *CCNzy,
                           cufftComplex *CCNzz,
                           cufftComplex *CCHx,
                           cufftComplex *CCHy,
                           cufftComplex *CCHz,
                           cufftComplex *CCSx,
                           cufftComplex *CCSy,
                           cufftComplex *CCSz
                           )
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            CCHx[i].x = (CCNxx[i].x*CCSx[i].x - CCNxx[i].y*CCSx[i].y + CCNxy[i].x*CCSy[i].x - CCNxy[i].y*CCSy[i].y + CCNxz[i].x*CCSz[i].x - CCNxz[i].y*CCSz[i].y);
            CCHx[i].y = (CCNxx[i].x*CCSx[i].y + CCNxx[i].y*CCSx[i].x + CCNxy[i].x*CCSy[i].y + CCNxy[i].y*CCSy[i].x + CCNxz[i].x*CCSz[i].y + CCNxz[i].y*CCSz[i].x);

            CCHy[i].x = (CCNyx[i].x*CCSx[i].x - CCNyx[i].y*CCSx[i].y + CCNyy[i].x*CCSy[i].x - CCNyy[i].y*CCSy[i].y + CCNyz[i].x*CCSz[i].x - CCNyz[i].y*CCSz[i].y);
            CCHy[i].y = (CCNyx[i].x*CCSx[i].y + CCNyx[i].y*CCSx[i].x + CCNyy[i].x*CCSy[i].y + CCNyy[i].y*CCSy[i].x + CCNyz[i].x*CCSz[i].y + CCNyz[i].y*CCSz[i].x);

            CCHz[i].x = (CCNzx[i].x*CCSx[i].x - CCNzx[i].y*CCSx[i].y + CCNzy[i].x*CCSy[i].x - CCNzy[i].y*CCSy[i].y + CCNzz[i].x*CCSz[i].x - CCNzz[i].y*CCSz[i].y);
            CCHz[i].y = (CCNzx[i].x*CCSx[i].y + CCNzx[i].y*CCSx[i].x + CCNzy[i].x*CCSy[i].y + CCNzy[i].y*CCSy[i].x + CCNzz[i].x*CCSz[i].y + CCNzz[i].y*CCSz[i].x);
        }
    }

    __global__ void CBfDip(int N,
                           double MsV,
                           float *CHDemag,
                           float *Cfspin,
                           float *Ccoord)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            float ri[3]={Ccoord[3*i],Ccoord[3*i+1],Ccoord[3*i+2]};
            register float locField[3]={0,0,0};
            register unsigned int j=0;
            for(j = 0 ; j < N ; j++)
            {
                if(i!=j)
                {
                    float rj[3]={Ccoord[3*j],Ccoord[3*j+1],Ccoord[3*j+2]};
                    float rij[3]={ri[0]-rj[0],ri[1]-rj[1],ri[2]-rj[2]};
                    //float modR=sqrt( rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2] );
                    float ooSqrtR=rsqrtf(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                    //float ooR3=1./(modR*modR*modR);
                    float ooR3=ooSqrtR*ooSqrtR*ooSqrtR;
                    float eij[3]={rij[0]*ooSqrtR,rij[1]*ooSqrtR,rij[2]*ooSqrtR};

                    float fMsV=float(MsV);
                    float m[3]={Cfspin[3*j]*fMsV,Cfspin[3*j+1]*fMsV,Cfspin[3*j+2]*fMsV};
                    float m_dot_eij=m[0]*eij[0]+m[1]*eij[1]+m[2]*eij[2];

                    locField[0]-=(m[0]-3.0*m_dot_eij*eij[0])*ooR3;
                    locField[1]-=(m[1]-3.0*m_dot_eij*eij[1])*ooR3;
                    locField[2]-=(m[2]-3.0*m_dot_eij*eij[2])*ooR3;

                }
            }
            //apply the 4pi/mu0 prefactor and set the field in global memory
            for(j = 0 ; j < 3 ; j++)
            {
                locField[j]*=1e-7;
                CHDemag[3*i+j]=locField[j];
            }
        }
    }
    __global__ void CZeroCrand(int N,
                           float *Crand
                           )
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            for(unsigned int j = 0 ; j < 6 ; j++)
            {
                Crand[6*i+j]=0.0;
            }
        }
    }
    __global__ void CZeroCNSpinArray(int N,
                           float *Car
                           )
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            Car[i]=0.0;
        }
    }
    //This needs to be done with a seperate kernel because the size (N)
    //of the zero padded spin arrays is bigger than the number of spins
    __global__ void CCopySpin(int zpN,unsigned int N,double *Cspin,int *Czpsn,cufftComplex *CCSx,cufftComplex *CCSy,cufftComplex *CCSz)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<zpN)
        {
            CCSx[i].x=0;CCSx[i].y=0;
            CCSy[i].x=0;CCSy[i].y=0;
            CCSz[i].x=0;CCSz[i].y=0;
            if(i<N)
            {
                //lookup the array value for spin i in the zero pad array
                int lzpsn=Czpsn[i];
                //copy the spin data to the zero padded spin arrays
                //for the fourier transform
                CCSx[lzpsn].x=Cspin[3*i];
                CCSy[lzpsn].x=Cspin[3*i+1];
                CCSz[lzpsn].x=Cspin[3*i+2];
            }
        }
    }

    __global__ void CCopyFields(int N,int zpN,float *CHDemag,int *Czpsn,cufftComplex *CCHx,cufftComplex *CCHy,cufftComplex *CCHz)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            int lzpsn=Czpsn[i];
            CHDemag[3*i]=float(CCHx[lzpsn].x)/float(zpN);
            CHDemag[3*i+1]=float(CCHy[lzpsn].x)/float(zpN);
            CHDemag[3*i+2]=float(CCHz[lzpsn].x)/float(zpN);
        }
    }
}
