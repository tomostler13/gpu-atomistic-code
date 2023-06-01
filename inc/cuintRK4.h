#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#ifndef _CUINTRK4_H_
#define _CUINTRK4_H_
namespace cuintRK4
{
    extern __global__ void CdetRK4k1(int,double,double,double,double,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK4k2(int,double,double,double,double,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK4k3(int,double,double,double,double,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK4k4(int,double,double,double,double,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern void copyConstData();
    extern __global__ void CdetSnp1(int,double*,double*,double*,double*,double*);
}
#endif /*_CUINTRK4_H_*/
