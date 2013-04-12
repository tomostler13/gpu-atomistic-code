// File: cuda.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 12 Apr 2013 14:43:41
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <stdio.h>
#include "../inc/arrays.h"
#ifndef _CULLG_H_
#define _CULLG_H_
namespace cullg
{
    extern cudaDeviceProp deviceProp;
    extern curandGenerator_t gen;
    extern int nrank;
    extern bool initosT;
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
    void llgGPU(unsigned int&,Array<double>&);
}
#endif /*_CULLB_H_*/
