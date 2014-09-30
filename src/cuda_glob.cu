// File: cuda.cu
// Author:Tom Ostler
// Created: 26/06/2014
// Last-modified: 30 Sep 2014 19:31:54
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/intmat.h"
#include "../inc/util.h"
#include "../inc/fields.h"
#include "../inc/arrays.h"
#include "../inc/cudadefs.h"
#include "../inc/defines.h"
#include "../inc/cufields.h"
#include "../inc/cuint.h"
#include "../inc/llg.h"
//Cuda headers
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cufft.h>
//Library headers
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
namespace cullg
{
    cudaDeviceProp deviceProp;
    curandGenerator_t gen;
    //number of threads per block and blocks per grid
    int threadsperblock,blockspergrid;
    //same but for the zero padded work spaces
    int zpblockspergrid;
    //same but for the complex zero padded work space (N_Z/2+1) for r2c and c2r transforms
    int czpblockspergrid;
    //rank of the FFT
    int nrank=3;
    //device pointers for Fourier space calculations
    cufftComplex *CNk=NULL;
    cufftComplex *CSk=NULL;
    cufftReal *CSr=NULL;
    cufftComplex *CHk=NULL;
    cufftReal *CHr=NULL;
    //unsigned int the kx, ky and kz positions of the spins. The point is that you can use these arrays to
    //lookup which element of the array the the spin data should be copied to.
    unsigned int *Ckx=NULL,*Cky=NULL,*Ckz=NULL,*Cspec=NULL;

    //device pointers
    double *Cspin=NULL,*Cespin=NULL,*Crand=NULL;
    float *CH=NULL;
    double *Clambda=NULL,*Csigma=NULL,*Cfn=NULL,*Cllgpf=NULL;
    //cufft plans
    cufftHandle C3DPr2c,C3DPc2r;

}
