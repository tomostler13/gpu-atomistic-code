// File: cuda.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 22 Nov 2014 16:58:32
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <stdio.h>
#ifndef _CULLG_H_
#define _CULLG_H_
namespace cullg
{
    extern cudaDeviceProp deviceProp;
    extern curandGenerator_t gen;
    //number of threads per block and blocks per grid
    extern int threadsperblock,blockspergrid;
    //same but for the zero padded work spaces
    extern int zpblockspergrid;
    //This is the number of block per grid for addressing the elements of the real space
    //spin and field arrays (dimensions: NUMSPEC x 3 x ZPDIM[0] x ZPDIM[1] x ZPDIM[2]
    extern int rsarzpblockspergrid;
    //rank of the FFT
    extern int nrank;
    //lookup for the x,y and z locations (on k-mesh) of the spins
    extern unsigned int *Ckx,*Cky,*Ckz,*Cspec,*Cxadj,*Cadjncy;
    //device pointers for Fourier space calculations
    extern cufftComplex *CNk,*CSk,*CSr,*CHk,*CHr;
    //offsets for the DIA sparse matrix multiplication
    extern int *Cdiagoffset,*Coffdiagoffset;
    //exchange tensor in DIA format
    extern float *Cdxx,*Cdyy,*Cdzz;
    extern float *Cdxy,*Cdxz,*Cdyx,*Cdyz,*Cdzx,*Cdzy;

    //device pointers
    extern  float *Crand,*CH,*CHDemag,*Cmagmom;
    extern  double *Cfn,*Csigma,*Clambda,*Cllgpf,*Cspin,*Cespin,*Ck1u,*Ck1udir;
    //cufft plans
    extern cufftHandle SPc2c,FPc2c;
    //device pointers
    void allocate_memory_on_card();
    void setup_fourier_transform();
    void spins_forward();
    void fields_back();
    void initGPU();
    inline void check_cuda_errors(const char *filename, const int line_number)
    {
        #ifdef DEBUG
        cudaThreadSynchronize();
        cudaError_t error = cudaGetLastError();
        if(error != cudaSuccess)
        {
            printf("CUDA error at %s:%i: %s\n", filename, line_number, cudaGetErrorString(error));
            exit(-1);
        }
        #endif
    };
    void cuinit(int argc,char *argv[]);
    void deallocate_cuda_memory();
    void llgGPU(unsigned int&);
}
#endif /*_CULLB_H_*/
