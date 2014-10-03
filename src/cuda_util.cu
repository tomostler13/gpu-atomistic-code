// File: cuda.cu
// Author:Tom Ostler
// Created: 26/06/2014
// Last-modified: 03 Oct 2014 14:43:49
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
//The function of this file is to house a number of routines
//that deal with a number of underlying routines, such as
//mallocing/de(m)allocing memory, setting up fft's etc.
// Requires: cullg::cuinit() to be called
namespace cullg
{

    void setup_fourier_transform()
    {
        /*Create a 3D FFT plan. */
        int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
        int *inembed=n;
        int *onembed=n;
        int istride=1;
        int ostride=1;
        int idist=geom::zps;
        int odist=geom::zps;
        config::openLogFile();
        config::printline(config::Log);
        FIXOUT(config::Log,"Parameters entering into CUFFT plan of the spin arrays (forward)" << std::endl);
        FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
        FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
        FIXOUT(config::Log,"How many (FFT's) = " << geom::ucm.GetNMS()*3 << std::endl);
        FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
        FIXOUT(config::Log,"istride = " << istride << std::endl);
        FIXOUT(config::Log,"idist = " << idist << std::endl);
        FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
        FIXOUT(config::Log,"ostride = " << ostride << std::endl);
        FIXOUT(config::Log,"odist = " << odist << std::endl);
        FIXOUT(config::Log,"Direction (sign) = " << "CUFFTW_FORWARD" << std::endl);
        if(cufftPlanMany(&SPc2c,3,n,inembed,istride,idist,onembed,ostride,odist,CUFFT_C2C,geom::ucm.GetNMS()*3)!=CUFFT_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CUFFT 3D plan creation failed");
        }
        if(cufftPlanMany(&FPc2c,3,n,onembed,ostride,odist,inembed,istride,idist,CUFFT_C2C,geom::ucm.GetNMS()*3)!=CUFFT_SUCCESS)
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
                                for(unsigned int k = 0 ; k < geom::zpdim[1]*geom::Nk[2] ; k++)
                                {
                                    for(unsigned int l = 0 ; l < 2 ; l++)
                                    {
                                        tempNkab(s1,s2,alpha,beta,i,j,k)[l]=static_cast<float>(intmat::Nkab(s1,s2,alpha,beta,i,j,k)[l]);
                                    }

                                }
                            }
                        }
                    }
                }
            }
        }

        //copy the FT'd interaction matrix to the card
        CUDA_CALL(cudaMemcpy(CNk,tempNkab.ptr(),geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3*geom::zpdim[0]*geom::zpdim[1]*geom::zpdim[2]*geom::Nk[0]*geom::Nk[1]*geom::Nk[2]*sizeof(fftwf_complex),cudaMemcpyHostToDevice));
//        intmat::Nkab.clear();
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
        CUDA_CALL(cudaFree(CSr));
        CUDA_CALL(cudaFree(CSk));
        CUDA_CALL(cudaFree(CHk));
        CUDA_CALL(cudaFree(CHr));
        CUDA_CALL(cudaFree(Cspin));
        CUDA_CALL(cudaFree(Cespin));
        CUDA_CALL(cudaFree(Crand));
        CUDA_CALL(cudaFree(CH));
        CUDA_CALL(cudaFree(Cfn));
        CUDA_CALL(cudaFree(Csigma));
        CUDA_CALL(cudaFree(Clambda));
        CUDA_CALL(cudaFree(Cllgpf));
        CUDA_CALL(cudaFree(Cspec));
        CUDA_CALL(cudaFree(Ckx));
        CUDA_CALL(cudaFree(Cky));
        CUDA_CALL(cudaFree(Ckz));
        config::Info << "Done" << std::endl;
    }
    void spins_forward()
    {
        Array5D<fftwf_complex> temp,temp1;
        temp.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        temp.IFill(0);
        temp1.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        temp1.IFill(0);
        CUDA_CALL(cudaMemcpy(temp.ptr(),CSr,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex),cudaMemcpyDeviceToHost));
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[0]*geom::Nk[0] ; k++)
                {
                    for(unsigned int l = 0 ; l < geom::zpdim[1]*geom::Nk[1] ; l++)
                    {
                        for(unsigned int m = 0 ; m < geom::zpdim[2]*geom::Nk[2] ; m++)
                        {
                            //std::cout << "rs-spins on CPU:\t" << temp(i,j,k,l,m)[0] << "\t" << temp(i,j,k,l,m)[1] << std::endl;

                        }
                    }
                }
            }
        }
        int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
        int *inembed=n;
        int *onembed=n;
        int istride=1;
        int ostride=1;
        int odist=geom::zps;
        int idist=geom::zps;
        fftwf_plan plan=fftwf_plan_many_dft(3,n,geom::ucm.GetNMS()*3,temp.ptr(),inembed,istride,idist,temp.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_PATIENT);
        fftwf_execute(plan);
        CUFFT_CALL(cufftExecC2C(SPc2c,CSr,CSk,CUFFT_FORWARD));
        CUDA_CALL(cudaMemcpy(temp1.ptr(),CSk,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex),cudaMemcpyDeviceToHost));
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[0]*geom::Nk[0] ; k++)
                {
                    for(unsigned int l = 0 ; l < geom::zpdim[1]*geom::Nk[1] ; l++)
                    {
                        for(unsigned int m = 0 ; m < geom::zpdim[2]*geom::Nk[2] ; m++)
                        {
                            if(fabs(temp(i,j,k,l,m)[0]-temp1(i,j,k,l,m)[0])>1e-4 || fabs(temp(i,j,k,l,m)[1]-temp1(i,j,k,l,m)[1])>1e-4)
                            {
                            //std::cout << "diff on CPU:\t" << temp(i,j,k,l,m)[0]-temp1(i,j,k,l,m)[0] << "\t" << temp(i,j,k,l,m)[1]-temp1(i,j,k,l,m)[1] << std::endl;
                            }

                        }
                    }
                }
            }
        }
        Array7D<fftwf_complex> temp2;
        temp2.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),3,3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        CUDA_CALL(cudaMemcpy(temp2.ptr(),CNk,geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3*geom::zps*sizeof(cufftComplex),cudaMemcpyDeviceToHost));
        temp.IFill(0);
        register unsigned int i = 0,j = 0, k = 0, s1 = 0, s2=0, alpha=0, beta=0;
        for(i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    for(s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
                    {
                        for(s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
                        {
                            for(alpha = 0 ; alpha < 3 ; alpha++)
                            {
                                for(beta = 0 ; beta < 3 ; beta++)
                                {
                                    temp(s1,alpha,i,j,k)[0]+=(temp2(s1,s2,alpha,beta,i,j,k)[0]*temp1(s2,beta,i,j,k)[0]-temp2(s1,s2,alpha,beta,i,j,k)[1]*temp1(s2,beta,i,j,k)[1]);
                                    temp(s1,alpha,i,j,k)[1]+=(temp2(s1,s2,alpha,beta,i,j,k)[0]*temp1(s2,beta,i,j,k)[1]+temp2(s1,s2,alpha,beta,i,j,k)[1]*temp1(s2,beta,i,j,k)[0]);
                                    if(fabs(temp2(s1,s2,alpha,beta,i,j,k)[0]-intmat::Nkab(s1,s2,alpha,beta,i,j,k)[0])>0.1 || fabs(temp2(s1,s2,alpha,beta,i,j,k)[1]-intmat::Nkab(s1,s2,alpha,beta,i,j,k)[1])>0.1)
                                    {
                                    std::cout << "GPU/CPU diff\t" << temp2(s1,s2,alpha,beta,i,j,k)[0]-intmat::Nkab(s1,s2,alpha,beta,i,j,k)[0] << "\t" << temp2(s1,s2,alpha,beta,i,j,k)[1]-intmat::Nkab(s1,s2,alpha,beta,i,j,k)[1] <<std::endl;
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }
        for(i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    for(s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
                    {
                        for(s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
                        {
                            for(alpha = 0 ; alpha < 3 ; alpha++)
                            {
                                for(beta = 0 ; beta < 3 ; beta++)
                                {
                                    std::cout << temp(s1,alpha,i,j,k)[0] << "\t" << temp(s1,alpha,i,j,k)[1] << std::endl;
                                }
                            }
                        }
                    }

                }
            }
        }
        fftwf_plan plan1=fftwf_plan_many_dft(3,n,geom::ucm.GetNMS()*3,temp.ptr(),inembed,istride,idist,temp.ptr(),onembed,ostride,odist,FFTW_BACKWARD,FFTW_PATIENT);
        fftwf_execute(plan1);
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[0]*geom::Nk[0] ; k++)
                {
                    for(unsigned int l = 0 ; l < geom::zpdim[1]*geom::Nk[1] ; l++)
                    {
                        for(unsigned int m = 0 ; m < geom::zpdim[2]*geom::Nk[2] ; m++)
                        {
                            std::cout << "k-fields on CPU:\t" << temp(i,j,k,l,m)[0] << "\t" << temp(i,j,k,l,m)[1] << std::endl;

                        }
                    }
                }
            }
        }


    }

    void fields_back()
    {
        Array5D<fftwf_complex> temp;
        temp.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        temp.IFill(0);
        CUDA_CALL(cudaMemcpy(temp.ptr(),CHk,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex),cudaMemcpyDeviceToHost));
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[0]*geom::Nk[0] ; k++)
                {
                    for(unsigned int l = 0 ; l < geom::zpdim[1]*geom::Nk[1] ; l++)
                    {
                        for(unsigned int m = 0 ; m < geom::zpdim[2]*geom::Nk[2] ; m++)
                        {
                            std::cout << "k-fields on CPU:\t" << temp(i,j,k,l,m)[0] << "\t" << temp(i,j,k,l,m)[1] << std::endl;

                        }
                    }
                }
            }
        }
        CUFFT_CALL(cufftExecC2C(FPc2c,CHk,CHr,CUFFT_INVERSE));
    }

    void allocate_memory_on_card()
    {
        //all of the GPU memory allocations should happen here.
        //--------------------------------------------------------------------------------
        CUDA_CALL(cudaMalloc((void**)&CNk,geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3*geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CSk,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CSr,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CHk,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CHr,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&Cspin,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Cespin,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Crand,3*geom::nspins*sizeof(float)));
        CUDA_CALL(cudaMalloc((void**)&CH,3*geom::nspins*sizeof(float)));
        CUDA_CALL(cudaMalloc((void**)&Cfn,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Clambda,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Csigma,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Cllgpf,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Ckx,geom::nspins*sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc((void**)&Cky,geom::nspins*sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc((void**)&Ckz,geom::nspins*sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc((void**)&Cspec,geom::nspins*sizeof(unsigned int)));
        //--------------------------------------------------------------------------------
        //this section sorts out the copying of the data from the CPU to the card
        //--------------------------------------------------------------------------------
        //copy the sigma prefactor
        CUDA_CALL(cudaMemcpy(Csigma,geom::sigma.ptr(),geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Clambda,geom::lambda.ptr(),geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Cllgpf,geom::llgpf.ptr(),geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        //declare some arrays for doing copying to card
        //Nspins float array, 3*Nspins float array.
        float *nsfa=new float[geom::nspins];
        float *tnsfa=new float[3*geom::nspins];
        //Nspins double array, 3*Nspins double array
        double *nsda=new double[geom::nspins];
        double *tnsda=new double[3*geom::nspins];
        //Nspins int array, 3*Nspins int array
        int *nsia=new int[geom::nspins];
        //copy the location of the spins in real space to the device
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            nsia[i]=geom::lu(i,0);
        }
        CUDA_CALL(cudaMemcpy(Ckx,nsia,geom::nspins*sizeof(unsigned int),cudaMemcpyHostToDevice));
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            nsia[i]=geom::lu(i,1);
        }
        CUDA_CALL(cudaMemcpy(Cky,nsia,geom::nspins*sizeof(unsigned int),cudaMemcpyHostToDevice));
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            nsia[i]=geom::lu(i,2);
        }
        CUDA_CALL(cudaMemcpy(Ckz,nsia,geom::nspins*sizeof(unsigned int),cudaMemcpyHostToDevice));
        //and copy the species list
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            nsia[i]=geom::lu(i,3);
        }
        CUDA_CALL(cudaMemcpy(Cspec,nsia,geom::nspins*sizeof(unsigned int),cudaMemcpyHostToDevice));
        int *tnsia=new int[3*geom::nspins];
        //copy spin data to single array
        util::copy3vecto1(geom::nspins,spins::Sx,spins::Sy,spins::Sz,tnsda);
        //copy spin data to card
        CUDA_CALL(cudaMemcpy(Cspin,tnsda,3*geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        //zero the field array
        for(unsigned int i = 0 ; i < 3*geom::nspins ; i++){tnsfa[i]=0.0;}CUDA_CALL(cudaMemcpy(CH,tnsfa,3*geom::nspins*sizeof(float),cudaMemcpyHostToDevice));
        //call the kernel to zero the spin array


        //make sure we clean up when the program exits
        atexit(deallocate_cuda_memory);
    }
}
