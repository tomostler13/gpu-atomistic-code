#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#ifndef _CUINTRK5_H_
#define _CUINTRK5_H_
namespace cuintRK5
{
    extern __global__ void CdetRK5k1(int,double,double,double,double,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK5k2(int,double,double,double,double,double*,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK5k3(int,double,double,double,double,double*,double*,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK5k4(int,double,double,double,double,double*,double*,double*,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK5k5(int,double,double,double,double,double*,double*,double*,double*,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern __global__ void CdetRK5k6(int,double,double,double,double,double*,double*,double*,double*,double*,double*,double*,double*,double*,float*,double*,double*,double*,double*,double*,double*);
    extern void copyConstData();
    extern __global__ void CdetSnp1(int,double*,double*,double*,double*,double*,double*,double*);
}
#endif /*_CUINTRK5_H_*/
