#include <cuda.h>
#include <cufft.h>
#ifndef _CUFIELDS_H_
#define _CUFIELDS_H_
namespace cufields
{
    extern __global__ void CFConv(int,
                                  cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,
                                  cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,cufftComplex*,
                                  cufftComplex*,cufftComplex*,cufftComplex*);
    extern __global__ void CCopySpin(int,unsigned int,double*,int*,cufftReal*,cufftReal*,cufftReal*,cufftReal*,cufftReal*,cufftReal*);
    extern __global__ void CCopyFields(int,int,float*,int*,cufftReal*,cufftReal*,cufftReal*);
    extern __global__ void CBFDip(int,float*,double*,double*);
    extern __global__ void CZeroField(int,float*);
}
#endif /*_CUFIELDS_H_*/
