#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#ifndef _CUINT_H_
#define _CUINT_H_
#define NR_END 1
#define FREE_ARG char*
#define FREERETURN {free_dmatrix(fjac,1,NLAT,1,NLAT);free_dvector(fvec,1,NLAT);\
	free_dvector(p,1,NLAT);free_ivector(indx,1,NLAT);return;}

namespace cuint
{
    extern __global__ void CSetFunctionPointer(int,int,int,int,int,int,int,int,unsigned int,unsigned int);
    extern __global__ void CHeun1(int,float,float,float,double,float*,float *,double *,double *,float *,int *,int *,double *,double *,float *,double *);
    extern __global__ void CHeun2(int,float,float,float,double,float*,double*,float *,float *,float *,double *,double *,int *,int *,double *,double *);
    void copyConstData();


}
#endif /*_CUINT_H_*/
