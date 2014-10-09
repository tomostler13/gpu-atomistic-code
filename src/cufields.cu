// File: cufields.cu
// Author:Tom Ostler
// Last-modified: 09 Oct 2014 14:29:15
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
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/llg.h"
namespace cufields
{
    //constant memory declarations
    //The number of k-points and the zero pad size
    __constant__ unsigned int K[3]={0,0,0},ZPDIM[3]={0,0,0};
    //DND = Diagonal number of off diagonals, ODND = off-diagonal number of off diagonals
    __constant__ unsigned int DND=0,ODND=0;

    void copyConstData()
    {
        FIXOUT(config::Info,"Copying const data with cufields scope to card:" << std::flush);
        cudaMemcpyToSymbol(K,geom::Nk.ptr(),3*sizeof(unsigned int));
        cudaMemcpyToSymbol(ZPDIM,&geom::zpdim,3*sizeof(unsigned int));
        cudaMemcpyToSymbol(DND,&exch::diagnumdiag,sizeof(unsigned int));
        cudaMemcpyToSymbol(ODND,&exch::offdiagnumdiag,sizeof(unsigned int));
        config::Info << "Done" << std::endl;
    }


    //perform the DIA matrix multiplication. Algorithm taken from:
    // Efficient Spart Matrix-Vector Multiplication on CUDA
    // Nathan Bell and Michael Garland
    // http://www.nvidia.com/docs/IO/66889/nvr-2008-004.pdf
    __global__ void CSpVM_DIA(int N,
            int *offset,
            float *dataxx,float *datayy,float *datazz,
            double *Cspin,
            float *CHDemag,float *CH)
    {
        const int row = blockDim.x*blockIdx.x + threadIdx.x;
        //num_rows=N as we are ALWAYS dealing with a num_rows=num_cols
        if(row < N)
        {
            //contains sum over neighbours
            float dot[3]={0,0,0};
            for(int n = 0 ; n < DND ; n++)
            {
                int col = row + offset[n];
                float val[3] = {dataxx[N*n + row],datayy[N*n + row],datazz[N*n + row]};
                if(col >= 0 && col < N)
                {
                    for(unsigned int co = 0 ; co < 3 ; co++)
                    {
                        dot[co]+=(val[co]*Cspin[3*col+co]);
                    }
                }
            }
            CH[3*row]=CHDemag[3*row]+dot[0];
            CH[3*row+1]=CHDemag[3*row+1]+dot[1];
            CH[3*row+2]=CHDemag[3*row+2]+dot[2];

        }
    }
    //perform the convolution in Fourier space
    __global__ void CFConv(int N,
            unsigned int NMS,
            cufftComplex *CNk,
            cufftComplex *CHk,
            cufftComplex *CSk
            )
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //the number of threads is the zps (zero pad size). We can then find the coordinate of the
            //fourier space k-point
            const unsigned int kx=i/(ZPDIM[1]*K[1]*ZPDIM[2]*K[2]),ky=i%(ZPDIM[1]*K[1]*ZPDIM[2]*K[2])/(ZPDIM[2]*K[2]),kz=i%(ZPDIM[2]*K[2]);

            for(unsigned int s1 = 0 ; s1 < NMS ; s1++)
            {
                for(unsigned int s2 = 0 ; s2 < NMS ; s2++)
                {
                    for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                    {

                        unsigned int hfari=(((s1*3+alpha)*ZPDIM[0]*K[0]+kx)*ZPDIM[1]*K[1]+ky)*ZPDIM[2]*K[2]+kz;
                        for(unsigned int beta = 0 ; beta < 3 ; beta++)
                        {
                            //calculate the interaction matrix array element
                            //from the 7D array lookup
                            //(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o
                            //dim0 = NMS, dim1 = NMS, dim2 = 3, dim3 = 3,
                            //dim4 = ZPDIM[0], dim5 = ZPDIM[1], dim6 = ZPDIM[2]

                            unsigned int Nari=(((((s1*NMS+s2)*3+alpha)*3+beta)*ZPDIM[0]*K[0]+kx)*ZPDIM[1]*K[1]+ky)*ZPDIM[2]*K[2]+kz;

                            //Calculate the field and spin array element (5D lookup)
                            //(((i*dim1+j)*dim2+k)*dim3+l)*dim4+m
                            //dim0 = NMS , dim1 = 3
                            //dim2 = ZPDIM[0], dim3 = ZPDIM[1], dim4 = ZPDIM[2]
                            unsigned int sfari=(((s2*3+beta)*ZPDIM[0]*K[0]+kx)*ZPDIM[1]*K[1]+ky)*ZPDIM[2]*K[2]+kz;
                            CHk[hfari].x += (CNk[Nari].x*CSk[sfari].x - CNk[Nari].y*CSk[sfari].y);
                            CHk[hfari].y += (CNk[Nari].x*CSk[sfari].y + CNk[Nari].y*CSk[sfari].x);
                        }
                    }
                }
            }
        }
    }
    //perform the convolution in Fourier space
    __global__ void CdipFConv(int N,
                           cufftComplex *CNk,
                           cufftComplex *CHk,
                           cufftComplex *CSk
                           )
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //the number of threads is the zps (zero pad size). We can then find the coordinate of the
            //fourier space k-point
            const unsigned int kx=i/(ZPDIM[1]*K[1]*ZPDIM[2]*K[2]),ky=i%(ZPDIM[1]*K[1]*ZPDIM[2]*K[2])/(ZPDIM[2]*K[2]),kz=i%(ZPDIM[2]*K[2]);

            for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
            {

                unsigned int hfari=((alpha*ZPDIM[0]*K[0]+kx)*ZPDIM[1]*K[1]+ky)*ZPDIM[2]*K[2]+kz;
                for(unsigned int beta = 0 ; beta < 3 ; beta++)
                {
                    //calculate the interaction matrix array element
                    //from the 5D array lookup
                    //(((i*dim1+j)*dim2+k)*dim3+l)*dim4+m
                    //dim0 = 3, dim1 = 3, dim2 = ZPDIM[0]*K[0], dim3 = ZPDIM[1]*K[1],
                    //dim4 = ZPDIM[2]*K[2]

                    unsigned int Nari=(((alpha*3+beta)*ZPDIM[0]*K[0]+kx)*ZPDIM[1]*K[1]+ky)*ZPDIM[2]*K[2]+kz;

                    //Calculate the field and spin array element (5D lookup)
                    //((i*dim1+j)*dim2+k)*dim3+l
                    //dim0 = 3 , dim1 = ZPDIM[0], dim2 = ZPDIM[1], dim3 = ZPDIM[2]
                    unsigned int sfari=((beta*ZPDIM[0]*K[0]+kx)*ZPDIM[1]*K[1]+ky)*ZPDIM[2]*K[2]+kz;
                    CHk[hfari].x += (CNk[Nari].x*CSk[sfari].x - CNk[Nari].y*CSk[sfari].y);
                    CHk[hfari].y += (CNk[Nari].x*CSk[sfari].y + CNk[Nari].y*CSk[sfari].x);
                }
            }
        }
    }

    //This needs to be done with a seperate kernel because the size (N)
    //of the zero padded spin arrays is bigger than the number of spins
    __global__ void CdipCopySpin(int N,double *Cspin,cufftComplex *CSr,unsigned int *Ckx,unsigned int *Cky,unsigned int *Ckz,float *Cmagmom)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;

        if(i<N)
        {
            //For a 4D array lookup we need i,j,k,m,n
            //((i*dim1+j)*dim2+k)*dim3+l
            //For the spin arrays the indices correspond to
            // i -> spin component
            // j -> x-coordinate
            // k -> y-coordinate
            // l -> z-coordinate
            // Here lookup i,j,k,m (call them li,lj,lk,ll)
            unsigned int lj=Ckx[i],lk=Cky[i],ll=Ckz[i];
            //the magnetic moment of spin i
            float magmom=Cmagmom[i];
            //loop over the 3 spin coordinates (j)
            for(unsigned int li = 0 ; li < 3 ; li++)
            {
                unsigned int arlu=((li*ZPDIM[0]*K[0]+lj)*ZPDIM[1]*K[1]+lk)*ZPDIM[2]*K[2]+ll;
                CSr[arlu].x=static_cast<float>(Cspin[3*i+li])*magmom;
            }
        }
    }
   //This needs to be done with a seperate kernel because the size (N)
    //of the zero padded spin arrays is bigger than the number of spins
    __global__ void CCopySpin(int N,double *Cspin,cufftComplex *CSr,unsigned int *Ckx,unsigned int *Cky,unsigned int *Ckz,unsigned int *Cspec)
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
                unsigned int arlu=(((li*3+lj)*ZPDIM[0]*K[0]+lk)*ZPDIM[1]*K[1]+ll)*ZPDIM[2]*K[2]+lm;
                CSr[arlu].x=static_cast<float>(Cspin[3*i+lj]);
            }
        }
    }

    __global__ void CCopyFields(int N,int zpN,float *CH,cufftComplex *CHr,unsigned int *Ckx,unsigned int *Cky,unsigned int *Ckz,unsigned int *Cspec)
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
                unsigned int arluv=(((li*3+lj)*ZPDIM[0]*K[0]+lk)*ZPDIM[1]*K[1]+ll)*ZPDIM[2]*K[2]+lm;
                CH[3*i+lj]=(CHr[arluv].x)/static_cast<float>(zpN);
            }
        }
    }
    __global__ void CdipCopyFields(int N,int zpN,float *CH,cufftComplex *CHr,unsigned int *Ckx,unsigned int *Cky,unsigned int *Ckz)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //For a 5D array lookup we need i,j,k,m,n
            //((i*dim1+j)*dim2+k)*dim3+l
            //For the spin arrays the indices correspond to
            // i -> spin component
            // j -> x-coordinate
            // k -> y-coordinate
            // l -> z-coordinate
            // Here lookup i,k,l,m (call them li,lk,ll,lm
            unsigned int lj=Ckx[i],lk=Cky[i],ll=Ckz[i];
            //loop over the 3 spin coordinates (j)
            for(unsigned int li = 0 ; li < 3 ; li++)
            {
                unsigned int arluv=((li*ZPDIM[0]*K[0]+lj)*ZPDIM[1]*K[1]+lk)*ZPDIM[2]*K[2]+ll;
                CH[3*i+li]=(CHr[arluv].x)/static_cast<float>(zpN);
            }
        }
    }
    //Cuda Set to Zero 5D Real Space Arrays
    __global__ void CZero5DRSArrays(int N,cufftComplex *CHr,cufftComplex *CSr,cufftComplex *CHk,cufftComplex *CSk)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            CHr[i].x=0.0;
            CHr[i].y=0.0;
            CSr[i].x=0.0;
            CSr[i].y=0.0;
            CHk[i].x=0.0;
            CHk[i].y=0.0;
            CSk[i].x=0.0;
            CSk[i].y=0.0;
        }
    }
    //Cuda Set to Zero 4D Real Space Arrays
    __global__ void CZero4DRSArrays(int N,cufftComplex *CHr,cufftComplex *CSr,cufftComplex *CHk,cufftComplex *CSk)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            CHr[i].x=0.0;
            CHr[i].y=0.0;
            CSr[i].x=0.0;
            CSr[i].y=0.0;
            CHk[i].x=0.0;
            CHk[i].y=0.0;
            CSk[i].x=0.0;
            CSk[i].y=0.0;
        }
    }
}
