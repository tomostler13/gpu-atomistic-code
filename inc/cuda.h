// File: cuda.h
// Author:Tom Ostler
// Last-modified: 07 Jan 2013 16:17:51
// Formally cuLLB.h
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <stdio.h>
#ifndef _CULLB_H_
#define _CULLB_H_
namespace cullb
{
    extern bool cullbi;
    extern bool cullbinit;
    extern cudaDeviceProp deviceProp;
    extern curandGenerator_t gen;
    extern int nrank;
    //device pointers
    void allocate_memory_on_card(unsigned int&);
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
    void LLBGPU(unsigned int&);
}
#endif /*_CULLB_H_*/
