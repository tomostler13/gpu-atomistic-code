// File: cufields.cu
// Author:Tom Ostler
// Last-modified: 02 Oct 2014 14:46:04
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
#include "../inc/llg.h"
namespace cufields
{
    //This stores the interaction matrix dimensions for the lookup and the number
    //of species
    __constant__ unsigned int IMDIMS[7]={0,0,3,3,0,0,0},NUMSPEC=0;
    //The number of k-points and the zero pad size
    __constant__ unsigned int K[3]={0,0,0},ZPDIM[3]={0,0,0},CPLXDIM=0;
    //Reduced timestep
    __constant__ double Crdt;

    void copyConstData()
    {
        FIXOUT(config::Info,"Copying const data with cufields scope to card:" << std::flush);
        cudaMemcpyToSymbol(*(&Crdt),&llg::rdt,sizeof(double));
        cudaMemcpyToSymbol(K,geom::Nk.ptr(),3*sizeof(unsigned int));
        cudaMemcpyToSymbol(ZPDIM,&geom::zpdim,3*sizeof(unsigned int));
        unsigned int nms=geom::ucm.GetNMS();
        cudaMemcpyToSymbol(*(&NUMSPEC),&nms,sizeof(unsigned int));
        cudaMemcpyToSymbol(*(&CPLXDIM),&geom::cplxdim,sizeof(unsigned int));
        config::Info << "Done" << std::endl;
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
            const unsigned int kx=i/(ZPDIM[0]*K[0]*ZPDIM[1]*K[1]),ky=i%(ZPDIM[0]*K[0]*ZPDIM[1]*K[1])/(ZPDIM[2]*K[2]),kz=i%(ZPDIM[2]*K[2]);
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

                CSr[(((li*3+lj)*ZPDIM[0]*K[0]+lk)*ZPDIM[1]*K[1]+ll)*ZPDIM[2]*K[2]+lm].x=float(Cspin[3*i+lj]);
                //This can probably be removed
                CSr[(((li*3+lj)*ZPDIM[0]*K[0]+lk)*ZPDIM[1]*K[1]+ll)*ZPDIM[2]*K[2]+lm].y=0.0;
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
                CH[3*i+lj]=(CHr[(((li*3+lj)*ZPDIM[0]*K[0]+lk)*ZPDIM[1]*K[1]+ll)*ZPDIM[2]*K[2]+lm].x)/static_cast<float>(zpN);
            }
        }
    }
    //Cuda Set to Zero 5D Real Space Arrays
    __global__ void CZero5DRSArrays(int N,cufftComplex *CHr,cufftComplex *CSr)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            CHr[i].x=0.0;
            CHr[i].y=0.0;
            CSr[i].x=0.0;
            CSr[i].y=0.0;
        }
    }
    //Cuda Set to Zero 5D Fourier Space Arrays
    __global__ void CZero5DFSArrays(int N,cufftComplex *CHk,cufftComplex *CSk)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            CHk[i].x=0.0;
            CHk[i].y=0.0;
            CSk[i].x=0.0;
            CSk[i].y=0.0;
        }
    }
}
