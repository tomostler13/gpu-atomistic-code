// File: cufields.cu
// Author:Tom Ostler
// Last-modified: 30 Sep 2014 19:58:53
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
    //This stores the interaction matrix dimensions for the lookup and the number
    //of species
    __constant__ unsigned int IMDIMS[7]={0,0,3,3,0,0,0},NUMSPEC=0;
    //The number of k-points and the zero pad size
    __constant__ unsigned int K[3]={0,0,0},ZPDIM[3]={0,0,0};
    //Reduced timestep
    __constant__ double Crdt;

    void copyConstData()
    {
        FIXOUT(config::Info,"Copying const data with cufields scope to card:" << std::flush);
        cudaMemcpyToSymbol(*(&Crdt),&llg::rdt,sizeof(double));
        cudaMemcpyToSymbol(K,geom::Nk.ptr(),3*sizeof(unsigned int));
        cudaMemcpyToSymbol(ZPDIM,&geom::zpdim,3*sizeof(unsigned int));
        cudaMemcpyToSymbol(*(&NUMSPEC),&geom::ucm.GetNMS(),sizeof(unsigned int));
        config::Info << "Done" << std::endl;
    }
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
    __global__ void CCopySpin(int N,double *Cspin,int *Czpsn,cufftReal *CSr,unsigned int *Ckx,unsigned int *Cky,unsigned int *Ckz,unsigned int *Cspec)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //For a 5D array lookup we need i,j,k,m,n
            //(((i*dim1+j)*dim2+k)*dim3+l)*dim4+m
            //For the spin arrays the indices correspond to
            // i -> species
            // j -> spin component
            // k -> x-coordinate
            // l -> y-coordinate
            // m -> z-coordinate
            // Here lookup i,k,l,m (call them li,lk,ll,lm
            unsigned int lk=Ckx[i],ll=Cky[i],lm=Ckz[i],li=Cspec[i];
            //loop over the 3 spin coordinates (j)
            for(unsigned int lj = 0 ; lj < 3 ; lj++)
            {
                unsigned int ac=
                CSr[(((li*3+lj)*ZPDIM[0]*K[0]+lk)*ZPDIM[1]*K[1]+ll)*ZPDIM[2]*K[2]+lm]=float(Cspin[3*i+lj]);
            }
        }
    }
}

    __global__ void CCopyFields(int N,int zpN,float *CH,int *Czpsn,cufftReal *CHr,unsigned int *Ckx,unsigned int *Cky,unsigned int *Ckz,unsigned int *Cspec)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //For a 5D array lookup we need i,j,k,m,n
            //(((i*dim1+j)*dim2+k)*dim3+l)*dim4+m
            //For the spin arrays the indices correspond to
            // i -> species
            // j -> spin component
            // k -> x-coordinate
            // l -> y-coordinate
            // m -> z-coordinate
            // Here lookup i,k,l,m (call them li,lk,ll,lm
            unsigned int lk=Ckx[i],ll=Cky[i],lm=Ckz[i],li=Cspec[i];
            //loop over the 3 spin coordinates (j)
            for(unsigned int lj = 0 ; lj < 3 ; lj++)
            {
                CH[3*i+lj]=(CCHx[(((li*3+lj)*ZPDIM[0]*K[0]+lk)*ZPDIM[1]*K[1]+ll)*ZPDIM[2]*K[2]+lm])/static_cast<float>(zpN);
            }
        }
    }
}
