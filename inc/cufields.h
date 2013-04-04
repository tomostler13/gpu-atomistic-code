#include <cuda.h>
#include <cufft.h>
#ifndef _CUFIELDS_H_
#define _CUFIELDS_H_
namespace cufields
{
    extern __global__ void CFConv(int,
                                  cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,
                                  cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*,
                                  cufftDoubleComplex*,cufftDoubleComplex*,cufftDoubleComplex*);
    extern __global__ void CCopySpin(int,unsigned int,double*,int*,cufftDoubleReal*,cufftDoubleReal*,cufftDoubleReal*,cufftDoubleReal*,cufftDoubleReal*,cufftDoubleReal*);
    extern __global__ void CCopyFields(int,int,double*,int*,cufftDoubleReal*,cufftDoubleReal*,cufftDoubleReal*);
}
#endif /*_CUFIELDS_H_*/
