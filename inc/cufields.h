#include <cuda.h>
#include <cufft.h>
#ifndef _CUFIELDS_H_
#define _CUFIELDS_H_
namespace cufields
{
    extern __global__ void CBfDip(int,double,float*,float*,float*);
    extern __global__ void CFConv(int,
                                  cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,
                                  cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,
                                  cufftComplex*,cufftComplex*,cufftComplex*);
    extern __global__ void CZeroCrand(int,float*);
    extern __global__ void CZeroNSpinArray(int,float*);
    extern __global__ void CCopySpin(int,unsigned int,double*,int*,cufftComplex*,cufftComplex*,cufftComplex*);
    extern __global__ void CCopyFields(int,int,float*,int*,cufftComplex*,cufftComplex*,cufftComplex*);
}
#endif /*_CUFIELDS_H_*/
