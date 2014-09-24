// File: cuda.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 23 Sep 2014 12:23:05
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
    //same but for the complex zero padded work space (N_Z/2+1) for r2c and c2r transforms
    extern int czpblockspergrid;
    //rank of the FFT
    extern int nrank;
    //device pointers for Fourier space calculations
    extern  cufftComplex *CCNxx;
    extern  cufftComplex *CCNxy;
    extern  cufftComplex *CCNxz;
    extern  cufftComplex *CCNyx;
    extern  cufftComplex *CCNyy;
    extern  cufftComplex *CCNyz;
    extern  cufftComplex *CCNzx;
    extern  cufftComplex *CCNzy;
    extern  cufftComplex *CCNzz;
    extern  cufftComplex *CCSkx;
    extern  cufftComplex *CCSky;
    extern  cufftComplex *CCSkz;
    extern  cufftReal *CCSrx;
    extern  cufftReal *CCSry;
    extern  cufftReal *CCSrz;
    extern  cufftComplex *CCHkx;
    extern  cufftComplex *CCHky;
    extern  cufftComplex *CCHkz;
    extern  cufftReal *CCHrx;
    extern  cufftReal *CCHry;
    extern  cufftReal *CCHrz;

    //device pointers
    extern  double *Cspin;
    extern  double *Cespin;
    extern  float *Crand;
    extern  float *CH;
    extern  unsigned int *Cspec;
    extern  int *Czpsn;//The is the zero pad spin number
    extern  int *Clu;
    extern  double *Cfn;
    //cufft plans
    extern cufftHandle C3DPr2c,C3DPc2r;
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
