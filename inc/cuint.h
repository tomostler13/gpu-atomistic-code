#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#ifndef _CUINT_H_
#define _CUINT_H_
namespace cuint
{
    //interaction matrix with uniform temperature
    extern __global__ void CHeun1(int,double,double*,double*,double*,double,float*,float*,double*,double*,float*,double*);
    //CSR neighbourlist with uniform temperature
    extern __global__ void CHeun1(int,double,double*,double*,double*,double,float*,float*,double*,double*,float*,double*,int*,int*,float*);
    //CSR neighbourlist with on-site temperature
    extern __global__ void CHeun1(int,double*,double*,double*,double*,double,float*,float*,double*,double*,float*,double*,int*,int*,float*);
    //interaction matrix with on-site temperature
    extern __global__ void CHeun1(int,double*,double*,double*,double*,double,float*,float*,double*,double*,float*,double*);
    //interaction matrix with uniform temperature
    extern __global__ void CHeun2(int,double,double*,double*,double*,double,float*,float*,double*,double*,float*,double*);
    //interaction matrix with on-site temperature
    extern __global__ void CHeun2(int,double*,double*,double*,double*,double,float*,float*,double*,double*,float*,double*);
    //CSR neighbourlist with uniform temperature
    extern __global__ void CHeun2(int,double,double*,double*,double*,double,float*,float*,double*,double*,float*,double*,int*,int*,float*);
    //CSR neighbourlist with on-site temperature
    extern __global__ void CHeun2(int,double*,double*,double*,double*,double,float*,float*,double*,double*,float*,double*,int*,int*,float*);
}
#endif /*_CUINT_H_*/
