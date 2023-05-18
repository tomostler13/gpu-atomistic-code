#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#ifndef _CUINT_H_
#define _CUINT_H_
namespace cuint
{
    extern __global__ void CHeun1(int,double,double,double,double,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CHeun2(int,double,double,double,double,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*,double*,double*);
    extern void copyConstData();
}
#endif /*_CUINT_H_*/
