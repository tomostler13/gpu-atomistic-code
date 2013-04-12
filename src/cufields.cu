// File: cufields.cu
// Author:Tom Ostler
// Last-modified: 12 Apr 2013 16:57:03
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
#include "../inc/cudadefs.h"
#include "../inc/defines.h"
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


    //This needs to be done with a seperate kernel because the size (N)
    //of the zero padded spin arrays is bigger than the number of spins
    __global__ void CCopySpin(int zpN,unsigned int N,double *Cspin,int *Czpsn,cufftReal *CCSx,cufftReal *CCSy,cufftReal *CCSz,cufftReal *CHrx,cufftReal *CHry,cufftReal *CHrz)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<zpN)
		{
			CCSx[i]=0.0;
			CCSy[i]=0.0;
			CCSz[i]=0.0;
			CHrx[i]=0.0;
			CHry[i]=0.0;
			CHrz[i]=0.0;
			//lookup the array value for spin i in the zero pad array
			int lzpsn=Czpsn[i];
			//copy the spin data to the zero padded spin arrays
			//for the fourier transform
			if(lzpsn>=0)
			{
				CCSx[i]=float(Cspin[3*lzpsn]);
				CCSy[i]=float(Cspin[3*lzpsn+1]);
				CCSz[i]=float(Cspin[3*lzpsn+2]);
			}
		}
    }

    __global__ void CCopyFields(int N,int zpN,float *CH,int *Czpsn,cufftReal *CCHx,cufftReal *CCHy,cufftReal *CCHz)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            int lzpsn=Czpsn[i];
            CH[3*i]=(CCHx[lzpsn])/float(zpN);
            CH[3*i+1]=(CCHy[lzpsn])/float(zpN);
            CH[3*i+2]=(CCHz[lzpsn])/float(zpN);
        }
    }

    //brute force calculation of the dipolar field
    //Only works for one species at the moment
    __global__ void CBFDip(int N,float *CH,double *Cr,double *Cspin)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            float ri[3]={Cr[3*i],Cr[3*i+1],Cr[3*i+2]};
            register double locField[3]={0,0,0};
            register unsigned int j=0;
            for(j = 0 ; j < N ; j++)
            {
                float rj[3]={Cr[3*j],Cr[3*j+1],Cr[3*j+2]};
                float Sj[3]={Cspin[3*j],Cspin[3*j+1],Cspin[3*j+2]};
                float rij[3]={ri[0]-rj[0],ri[1]-rj[1],ri[2]-rj[2]};
                float ooSqrtR=rsqrtf(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                float ooR3=ooSqrtR*ooSqrtR*ooSqrtR;
                float eij[3]={rij[0]*ooSqrtR,rij[1]*ooSqrtR,rij[2]*ooSqrtR};
                float s_dot_eij=Sj[0]*eij[0]+Sj[1]*eij[1]+Sj[2]*eij[2];

                for(unsigned int k = 0 ; k < 3 ; k++)
                {
                    locField[k]-=(Sj[k]-3.0*s_dot_eij*eij[k])*ooR3;
                }
            }
            //apply the 4pi/mu0 prefactor and set the field in global memory
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                locField[j]*=1e-7;
                CH[3*i+j]=locField[j];
            }

        }
    }
    __global__ void CZeroField(int N,float *CH)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            CH[3*i]=0.0;
            CH[3*i+1]=0.0;
            CH[3*i+2]=0.0;
        }
    }
}
