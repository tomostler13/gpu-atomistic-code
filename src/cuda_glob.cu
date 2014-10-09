// File: cuda.cu
// Author:Tom Ostler
// Created: 26/06/2014
// Last-modified: 09 Oct 2014 12:01:20
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
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
    //This is the number of block per grid for addressing the elements of the real space
    //spin and field arrays (dimensions: NUMSPEC x 3 x ZPDIM[0] x ZPDIM[1] x ZPDIM[2]
    int rsarzpblockspergrid=0;
    //Same as the rsarzpblockspergrid but with the ZPDIM[2] dimension now replaced with the (ZPDIM[2]+1)/2
    //for the c2r/r2c transforms
    int ksarzpblockspergrid=0;
    //rank of the FFT
    int nrank=3;
    //device pointers for Fourier space calculations
    cufftComplex *CNk=NULL;
    cufftComplex *CSk=NULL;
    cufftComplex *CSr=NULL;
    cufftComplex *CHk=NULL;
    cufftComplex *CHr=NULL;
    //unsigned int the kx, ky and kz positions of the spins. The point is that you can use these arrays to
    //lookup which element of the array the the spin data should be copied to.
    unsigned int *Ckx=NULL,*Cky=NULL,*Ckz=NULL,*Cspec=NULL;
    int *Cdiagoffset=NULL,*Coffdiagoffset=NULL;
    //DIA format components of the exchange tensor
    float *Cdxx=NULL,*Cdyy=NULL,*Cdzz=NULL;
    float *Cdxy=NULL,*Cdxz=NULL,*Cdyx=NULL,*Cdyz=NULL,*Cdzx=NULL,*Cdzy=NULL;

    //device pointers
    double *Cspin=NULL,*Cespin=NULL;
    float *CH=NULL,*Crand=NULL,*Cmagmom=NULL,*CHDemag=NULL;
    double *Clambda=NULL,*Csigma=NULL,*Cfn=NULL,*Cllgpf=NULL,*Ck1u=NULL,*Ck1udir=NULL;;
    //cufft plans
    cufftHandle FPc2c,SPc2c;

}
