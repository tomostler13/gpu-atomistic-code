// File: cuda.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 19 May 2023 13:24:15
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
    extern unsigned int *Ckx,*Cky,*Ckz,*Cspec,*Cxadj,*Cadjncy,*Cxadj_jkl,*Cadjncy_j,*Cadjncy_k,*Cadjncy_l;
    //device pointers for Fourier space calculations
    extern cufftComplex *CNk,*CSk,*CSr,*CHk,*CHr;
    //offsets for the DIA sparse matrix multiplication
    extern int *Cdiagoffset,*Coffdiagoffset;
    //exchange tensor in DIA format
    extern double *Cdxx,*Cdyy,*Cdzz;
    extern double *Cdxy,*Cdxz,*Cdyx,*Cdyz,*Cdzx,*Cdzy;

    //device pointers
    extern  float *Crand,*CHDemag,*Cmagmom;
    extern  double *CH,*Cfn,*Csigma,*Clambda,*Cllgpf,*Cspin,*Cespin,*Ck1u,*Ck1udir,*CDetFields,*CHstg,*CInitHstg;
    //for RK4 integration
    extern double *CRKk1,*CRKk2,*CRKk3,*CRKk4,*CRKk5,*CRKk6;
    //cufft plans
    extern cufftHandle SPc2c,FPc2c;
    extern int curandN;
    //device pointers
    void allocate_memory_on_card();
    void setup_fourier_transform();
    void spins_forward();
    void fields_back();
    void initGPU();
    inline void check_cuda_errors(const char *filename, const int line_number)
    {
        #ifdef DEBUG
        cudaDeviceSynchronize();
        cudaError_t error = cudaGetLastError();
        if(error != cudaSuccess)
        {
            printf("CUDA error at %s:%i: %s\n", filename, line_number, cudaGetErrorString(error));
            exit(-1);
        }
        #endif
    };
    void rampStag(double);
    void cuinit();
    void deallocate_cuda_memory();
    void CsetStagFieldsZero();
    void CsetStagFields();
    void llgGPU(unsigned int&);
    void llgGPURK4(unsigned int&);
    void llgGPURK5(unsigned int&);
}
#endif /*_CULLB_H_*/
