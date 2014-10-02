#include <cuda.h>
#include <cufft.h>
#ifndef _CUFIELDS_H_
#define _CUFIELDS_H_
namespace cufields
{
    extern __global__ void CFConv(int,unsigned int,cufftComplex*,cufftComplex*,cufftComplex*);
    extern __global__ void CCopySpin(int,double*,cufftReal*,unsigned int*,unsigned int*,unsigned int*,unsigned int*);
    extern __global__ void CCopyFields(int,int,float*,cufftReal*,unsigned int*,unsigned int*,unsigned int*,unsigned int*);
    extern __global__ void CZero5DRSArrays(int,cufftReal*,cufftReal*);
    extern __global__ void CZero5DFSArrays(int,cufftComplex*,cufftComplex*);
    extern void copyConstData();
}
#endif /*_CUFIELDS_H_*/
