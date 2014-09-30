// File: cuda.cu
// Author:Tom Ostler
// Created: 26/06/2014
// Last-modified: 30 Sep 2014 16:33:05
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
//The function of this file is to house a number of routines
//that deal with a number of underlying routines, such as
//mallocing/de(m)allocing memory, setting up fft's etc.
// Requires: cullg::cuinit() to be called
namespace cullg
{

    void setup_fourier_transform()
    {
        //Even though we have 9 interaction matrices, 3 field arrays and
        //3 spin arrays we only need one transform in cufft. This is because
        //we can reuse the plan and alternate the sign depending on whether
        //we have a forward or a back transform
        /*Create a 3D FFT plan. */
        if(cufftPlan3d(&C3DPr2c,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],CUFFT_R2C)!=CUFFT_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CUFFT 3D plan creation failed");
        }
        if(cufftPlan3d(&C3DPc2r,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],CUFFT_C2R)!=CUFFT_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CUFFT 3D plan creation failed");
        }


        //At this point we can copy the interaction matrix from the CPU
        //as there is no need to do the determination of the interaction
        //matrix on the card.
        //declare a holder on the heap
        Array7D<fftwf_complex> tempNkab;
        tempNkab.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),3,3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
        {
            for(unsigned int s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
            {
                for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                {
                    for(unsigned int beta = 0 ; beta < 3 ; beta++)
                    {
                        for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
                        {
                            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
                            {
                                for(unsigned int k = 0 ; k < geom::cplxdim ; k++)
                                {
                                    for(unsigned int l = 0 ; l < 2 ; l++)
                                    {
                                        tempNkab(s1,s2,alpha,beta,i,j,k)[l]=float(intmat::Nkab(s1,s2,alpha,beta,i,j,k)[l]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        //copy the FT'd interaction matrix to the card
        CUDA_CALL(cudaMemcpy(CNk,tempNkab.ptr(),geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3*geom::zpdim[0]*geom::zpdim[1]*geom::zpdim[2]*geom::Nk[0]*geom::Nk[1]*geom::Nk[2]*sizeof(cufftComplex),cudaMemcpyHostToDevice));
        intmat::Nkab.clear();
        //clear the floating point holding arrays as well
        tempNkab.clear();
        check_cuda_errors(__FILE__,__LINE__);
    }
    void deallocate_cuda_memory()
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**EXIT information***" << std::endl;
        FIXOUT(config::Info,"Freeing space on GPU device" << std::flush);
        CUDA_CALL(cudaFree(CNk));
        CUDA_CALL(cudaFree(CCSrx));
        CUDA_CALL(cudaFree(CCSry));
        CUDA_CALL(cudaFree(CCSrz));
        CUDA_CALL(cudaFree(CCSkx));
        CUDA_CALL(cudaFree(CCSky));
        CUDA_CALL(cudaFree(CCSkz));
        CUDA_CALL(cudaFree(CCHkx));
        CUDA_CALL(cudaFree(CCHky));
        CUDA_CALL(cudaFree(CCHkz));
        CUDA_CALL(cudaFree(CCHrx));
        CUDA_CALL(cudaFree(CCHry));
        CUDA_CALL(cudaFree(CCHrz));
        CUDA_CALL(cudaFree(Cspin));
        CUDA_CALL(cudaFree(Cespin));
        CUDA_CALL(cudaFree(Crand));
        CUDA_CALL(cudaFree(CH));
        CUDA_CALL(cudaFree(Czpsn));
        CUDA_CALL(cudaFree(Cfn));
        config::Info << "Done" << std::endl;
    }
    void spins_forward()
    {
        CUFFT_CALL(cufftExecR2C(C3DPr2c,CCSrx,CCSkx));
        CUFFT_CALL(cufftExecR2C(C3DPr2c,CCSry,CCSky));
        CUFFT_CALL(cufftExecR2C(C3DPr2c,CCSrz,CCSkz));
    }

    void fields_back()
    {
        CUFFT_CALL(cufftExecC2R(C3DPc2r,CCHkx,CCHrx));
        CUFFT_CALL(cufftExecC2R(C3DPc2r,CCHky,CCHry));
        CUFFT_CALL(cufftExecC2R(C3DPc2r,CCHkz,CCHrz));
    }

    void allocate_memory_on_card()
    {
        //all of the GPU memory allocations should happen here.
        //--------------------------------------------------------------------------------
        CUDA_CALL(cudaMalloc((void**)&CCNxx,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNxy,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNxz,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNyx,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNyy,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNyz,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNzx,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNzy,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNzz,geom::czps*sizeof(cufftComplex)));

        CUDA_CALL(cudaMalloc((void**)&CCSkx,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCSky,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCSkz,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCSrx,geom::zps*sizeof(cufftReal)));
        CUDA_CALL(cudaMalloc((void**)&CCSry,geom::zps*sizeof(cufftReal)));
        CUDA_CALL(cudaMalloc((void**)&CCSrz,geom::zps*sizeof(cufftReal)));
        CUDA_CALL(cudaMalloc((void**)&CCHkx,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCHky,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCHkz,geom::czps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCHrx,geom::zps*sizeof(cufftReal)));
        CUDA_CALL(cudaMalloc((void**)&CCHry,geom::zps*sizeof(cufftReal)));
        CUDA_CALL(cudaMalloc((void**)&CCHrz,geom::zps*sizeof(cufftReal)));
        CUDA_CALL(cudaMalloc((void**)&Cspin,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Cespin,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Crand,3*geom::nspins*sizeof(float)));
        CUDA_CALL(cudaMalloc((void**)&CH,3*geom::nspins*sizeof(float)));
        CUDA_CALL(cudaMalloc((void**)&Clu,geom::zps*sizeof(int)));
        CUDA_CALL(cudaMalloc((void**)&Czpsn,geom::nspins*sizeof(int)));
        CUDA_CALL(cudaMalloc((void**)&Cspec,geom::nspins*sizeof(unsigned int)));

        CUDA_CALL(cudaMalloc((void**)&Cfn,3*geom::nspins*sizeof(double)));

        //--------------------------------------------------------------------------------
        //this section sorts out the copying of the data from the CPU to the card
        //--------------------------------------------------------------------------------
        //declare some arrays for doing copying to card
        //Nspins float array, 3*Nspins float array.
        float *nsfa=new float[geom::nspins];
        float *tnsfa=new float[3*geom::nspins];
        //Nspins double array, 3*Nspins double array
        double *nsda=new double[geom::nspins];
        double *tnsda=new double[3*geom::nspins];
        //Nspins int array, 3*Nspins int array
        int *nsia=new int[geom::nspins];
        int *tnsia=new int[3*geom::nspins];
        //copy spin data to single array
        util::copy3vecto1(geom::nspins,spins::Sx,spins::Sy,spins::Sz,tnsda);
        //copy spin data to card
        CUDA_CALL(cudaMemcpy(Cspin,tnsda,3*geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        //copy the species information to the card
        CUDA_CALL(cudaMemcpy(Cspec,geom::spec,geom::nspins*sizeof(unsigned int),cudaMemcpyHostToDevice));
        int *sn=NULL;
        sn=new int[geom::zps];
        int count=0;
        for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    if(geom::coords(i,j,k,0) > -1)
                    {
                        sn[count]=geom::coords(i,j,k,0);
                    }
                    else
                    {
                        sn[count]=-1;
                    }
                    count++;
                }
            }
        }
        if(count!=geom::zps)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Error in counting zero pad size for GPU lookup");
        }
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            nsia[i]=spins::Srx.getarrayelement(geom::lu(i,0),geom::lu(i,1),geom::lu(i,2));
            //			std::cerr << nsia[i] << std::endl;
        }
        CUDA_CALL(cudaMemcpy(Czpsn,nsia,geom::nspins*sizeof(int),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Clu,sn,geom::zps*sizeof(int),cudaMemcpyHostToDevice));
        //zero the field array
        for(unsigned int i = 0 ; i < 3*geom::nspins ; i++){tnsfa[i]=0.0;}CUDA_CALL(cudaMemcpy(CH,tnsfa,3*geom::nspins*sizeof(float),cudaMemcpyHostToDevice));
        //		cufftReal
        //--------------------------------------------------------------------------------


        //make sure we clean up when the program exits
        atexit(deallocate_cuda_memory);
    }
}
