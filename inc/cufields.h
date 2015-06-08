#include <cuda.h>
#include <cufft.h>
#ifndef _CUFIELDS_H_
#define _CUFIELDS_H_
namespace cufields
{
    extern __global__ void CFConv(int,unsigned int,cufftComplex*,cufftComplex*,cufftComplex*);
    extern __global__ void CSpMV_DIA(int,int*,float*,float*,float*,double*,float*,float*);
    extern __global__ void CSpMV_CSR(unsigned int,unsigned int*,unsigned int*,float*,float*,float*,double*,float*,float*);
    extern __global__ void CdipFConv(int,cufftComplex*,cufftComplex*,cufftComplex*);
    extern __global__ void CCopySpin(int,double*,cufftComplex*,unsigned int*,unsigned int*,unsigned int*,unsigned int*);
    extern __global__ void CdipCopySpin(int,double*,cufftComplex*,unsigned int*,unsigned int*,unsigned int*,float*);
    extern __global__ void CCopyFields(int,int,float*,cufftComplex*,unsigned int*,unsigned int*,unsigned int*,unsigned int*);
    extern __global__ void CdipCopyFields(int,int,float*,cufftComplex*,unsigned int*,unsigned int*,unsigned int*);
    extern __global__ void CZero5DRSArrays(int,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*);
    extern __global__ void CZero4DRSArrays(int,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*);
    extern void copyConstData();
//    extern __global__ void CSpMV_CSR_FourSpin(unsigned int,unsigned int*,unsigned int*,unsigned int*,unsigned int*,float*,double*);
    extern __global__ void CSpMV_CSR_FourSpin(unsigned int,unsigned int*,unsigned int*,unsigned int*,unsigned int*,float*,double*);
}
#endif /*_CUFIELDS_H_*/
