// File: cuda.cu
// Author:Tom Ostler
// Created: 26/06/2014
// Last-modified: 18 May 2023 02:51:02 PM
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
#include "../inc/anis.h"
#include "../inc/exch.h"
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
        if(config::exchm==0 || (config::dipm==0 && config::inc_dip==true))
        {

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into CUFFT plan of the spin arrays (forward)" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
            int howmany=0;
            if(config::exchm==0)
            {
                howmany=geom::ucm.GetNMS()*3;
            }
            else if(config::exchm>0 && config::dipm==0 && config::inc_dip)
            {
                howmany=3;
            }
            FIXOUT(config::Log,"How many (FFT's) = " << howmany << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) = " << "CUFFT_FORWARD" << std::endl);
            if(cufftPlanMany(&SPc2c,3,n,inembed,istride,idist,onembed,ostride,odist,CUFFT_C2C,howmany)!=CUFFT_SUCCESS)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("CUFFT 3D plan creation failed");
            }
            else
            {
                FIXOUT(config::Log,"CUFFT returned success");
            }
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into CUFFT plan of the field arrays (inverse)" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << howmany << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",onembed[0],onembed[1],onembed[2]);
            FIXOUT(config::Log,"istride = " << ostride << std::endl);
            FIXOUT(config::Log,"idist = " << odist << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",inembed[0],inembed[1],inembed[2]);
            FIXOUT(config::Log,"ostride = " << istride << std::endl);
            FIXOUT(config::Log,"odist = " << idist << std::endl);
            FIXOUT(config::Log,"Direction (sign) = " << "CUFFT_INVERSE" << std::endl);
            if(cufftPlanMany(&FPc2c,3,n,onembed,ostride,odist,inembed,istride,idist,CUFFT_C2C,howmany)!=CUFFT_SUCCESS)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("CUFFT 3D plan creation failed");
            }
            else
            {
                FIXOUT(config::Log,"CUFFT returned success");
            }

            //At this point we can copy the interaction matrix from the CPU
            //as there is no need to do the determination of the interaction
            //matrix on the card.
            //declare a holder on the heap
            if(config::exchm==0)
            {
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
                                        for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
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
                CUDA_CALL(cudaMemcpy(CNk,tempNkab.ptr(),geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3*geom::zpdim[0]*geom::zpdim[1]*geom::zpdim[2]*geom::Nk[0]*geom::Nk[1]*geom::Nk[2]*sizeof(cufftComplex),cudaMemcpyHostToDevice));
                intmat::Nkab.clear();
                //clear the floating point holding arrays as well
                tempNkab.clear();
                check_cuda_errors(__FILE__,__LINE__);

            }
            else if(config::exchm>0 && config::dipm==0 && config::inc_dip==true)
            {
                Array5D<fftwf_complex> tempNkab;
                tempNkab.resize(3,3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
                tempNkab.IFill(0);
                for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                {
                    for(unsigned int beta = 0 ; beta < 3 ; beta++)
                    {
                        for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
                        {
                            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
                            {
                                for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                                {
                                    for(unsigned int l = 0 ; l < 2 ; l++)
                                    {
                                        tempNkab(alpha,beta,i,j,k)[l]=static_cast<float>(intmat::dipNkab(alpha,beta,i,j,k)[l]);
                                    }
                                }
                            }
                        }
                    }
                }
                //copy the FT'd interaction matrix to the card
                CUDA_CALL(cudaMemcpy(CNk,tempNkab.ptr(),3*3*geom::zpdim[0]*geom::zpdim[1]*geom::zpdim[2]*geom::Nk[0]*geom::Nk[1]*geom::Nk[2]*sizeof(cufftComplex),cudaMemcpyHostToDevice));
                //intmat::dipNkab.clear();
                //clear the floating point holding arrays as well
//                tempNkab.clear();
                check_cuda_errors(__FILE__,__LINE__);

            }
        }

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
        if(llg::intscheme==0)
        {
            CUDA_CALL(cudaFree(Cfn));
        }
        else if(llg::intscheme==1)
        {
            CUDA_CALL(cudaFree(CRK4k1));
            CUDA_CALL(cudaFree(CRK4k2));
            CUDA_CALL(cudaFree(CRK4k3));
            CUDA_CALL(cudaFree(CRK4k4));
        }
        CUDA_CALL(cudaFree(CDetFields));
        CUDA_CALL(cudaFree(Crand));
        CUDA_CALL(cudaFree(CH));
        CUDA_CALL(cudaFree(CHstg));
        CUDA_CALL(cudaFree(CInitHstg));
        CUDA_CALL(cudaFree(Csigma));
        CUDA_CALL(cudaFree(Clambda));
        CUDA_CALL(cudaFree(Cllgpf));
        CUDA_CALL(cudaFree(Cspec));
        CUDA_CALL(cudaFree(Ckx));
        CUDA_CALL(cudaFree(Cky));
        CUDA_CALL(cudaFree(Ckz));
        CUDA_CALL(cudaFree(CHDemag));
        CUDA_CALL(cudaFree(Cdxx));
        CUDA_CALL(cudaFree(Cdxy));
        CUDA_CALL(cudaFree(Cdxz));
        CUDA_CALL(cudaFree(Cdyx));
        CUDA_CALL(cudaFree(Cdyy));
        CUDA_CALL(cudaFree(Cdyz));
        CUDA_CALL(cudaFree(Cdzx));
        CUDA_CALL(cudaFree(Cdzy));
        CUDA_CALL(cudaFree(Cdzz));
        CUDA_CALL(cudaFree(Cmagmom));
        CUDA_CALL(cudaFree(Cdiagoffset));
        CUDA_CALL(cudaFree(Coffdiagoffset));
        config::Info << "Done" << std::endl;
    }
    void spins_forward()
    {
        CUFFT_CALL(cufftExecC2C(SPc2c,CSr,CSk,CUFFT_FORWARD));
    }

    void fields_back()
    {
        CUFFT_CALL(cufftExecC2C(FPc2c,CHk,CHr,CUFFT_INVERSE));
    }

    void allocate_memory_on_card()
    {
        //all of the GPU memory allocations should happen here.
        //--------------------------------------------------------------------------------
        if(config::exchm==0)
        {
            CUDA_CALL(cudaMalloc((void**)&CNk,geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3*geom::zps*sizeof(cufftComplex)));
            CUDA_CALL(cudaMalloc((void**)&CSk,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
            CUDA_CALL(cudaMalloc((void**)&CSr,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
            CUDA_CALL(cudaMalloc((void**)&CHk,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
            CUDA_CALL(cudaMalloc((void**)&CHr,geom::ucm.GetNMS()*3*geom::zps*sizeof(cufftComplex)));
            CUDA_CALL(cudaMalloc((void**)&Cspec,geom::nspins*sizeof(unsigned int)));
        }
        else if(config::exchm>0)
        {
            if(config::inc_dip==true && config::dipm==0)
            {
                CUDA_CALL(cudaMalloc((void**)&CNk,3*3*geom::zps*sizeof(cufftComplex)));
                CUDA_CALL(cudaMalloc((void**)&CSk,3*geom::zps*sizeof(cufftComplex)));
                CUDA_CALL(cudaMalloc((void**)&CSr,3*geom::zps*sizeof(cufftComplex)));
                CUDA_CALL(cudaMalloc((void**)&CHk,3*geom::zps*sizeof(cufftComplex)));
                CUDA_CALL(cudaMalloc((void**)&CHr,3*geom::zps*sizeof(cufftComplex)));
                CUDA_CALL(cudaMalloc((void**)&Cmagmom,geom::nspins*sizeof(float)));
            }

            if(config::exchm==1)//DIA
            {
                CUDA_CALL(cudaMalloc((void**)&Cdiagoffset,exch::diagoffset.size()*sizeof(int)));
                //copy the diag offsets
                CUDA_CALL(cudaMemcpy(Cdiagoffset,exch::diagoffset.ptr(),exch::diagoffset.size()*sizeof(int),cudaMemcpyHostToDevice));
                CUDA_CALL(cudaMalloc((void**)&Cdxx,exch::dataxx.size()*sizeof(double)));
                CUDA_CALL(cudaMalloc((void**)&Cdyy,exch::datayy.size()*sizeof(double)));
                CUDA_CALL(cudaMalloc((void**)&Cdzz,exch::datazz.size()*sizeof(double)));
                //dataxx, datayy and datazz should all be the same size
                double *temp=new double[exch::dataxx.size()];
                for(unsigned int i = 0 ; i < exch::dataxx.size() ; i++)
                {
                    temp[i]=exch::dataxx[i];
                }
                CUDA_CALL(cudaMemcpy(Cdxx,temp,exch::dataxx.size()*sizeof(double),cudaMemcpyHostToDevice));
                if(exch::dataxx.size()!=exch::datayy.size())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Something has gone wrong, the length of the diagonal (xx,yy and zz) components of the exchange\ntensor in the DIA format (exch::dataxx etc) should be the same. Here xx and yy not equal.");
                }
                for(unsigned int i = 0 ; i < exch::dataxx.size() ; i++)
                {
                    temp[i]=exch::datayy[i];
                }
                CUDA_CALL(cudaMemcpy(Cdyy,temp,exch::datayy.size()*sizeof(double),cudaMemcpyHostToDevice));
                if(exch::dataxx.size()!=exch::datazz.size())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Something has gone wrong, the length of the diagonal (xx,yy and zz) components of the exchange\ntensor in the DIA format (exch::dataxx etc) should be the same. Here xx and yy not equal.");
                }
                for(unsigned int i = 0 ; i < exch::dataxx.size() ; i++)
                {
                    temp[i]=exch::datayy[i];
                }
                CUDA_CALL(cudaMemcpy(Cdzz,temp,exch::datazz.size()*sizeof(double),cudaMemcpyHostToDevice));
                if(config::offdiag==true)
                {
                    CUDA_CALL(cudaMalloc((void**)&Coffdiagoffset,exch::offdiagoffset.size()*sizeof(int)));
                    CUDA_CALL(cudaMemcpy(Coffdiagoffset,exch::offdiagoffset.ptr(),exch::offdiagoffset.size()*sizeof(int),cudaMemcpyHostToDevice));
                    if(!(exch::dataxy.size()==exch::dataxz.size()==exch::datayx.size()==exch::datayz.size()==exch::datazx.size()==exch::datazy.size()))
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("Size of the off diagonal exchange components (in DIA format) are not of the same length.");
                    }
                    CUDA_CALL(cudaMalloc((void**)&Cdxy,exch::dataxy.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdxz,exch::dataxz.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdyx,exch::datayx.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdyz,exch::datayz.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdzx,exch::datazx.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdzy,exch::datazy.size()*sizeof(double)));
                    double *otemp=new double[exch::dataxy.size()];for(unsigned int i = 0 ; i < exch::dataxy.size() ; i++){otemp[i]=static_cast<double>(exch::dataxy[i]);}
                    CUDA_CALL(cudaMemcpy(Cdxy,otemp,exch::dataxy.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::dataxz.size() ; i++){otemp[i]=static_cast<double>(exch::dataxz[i]);}
                    CUDA_CALL(cudaMemcpy(Cdxz,otemp,exch::dataxz.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datayx.size() ; i++){otemp[i]=static_cast<double>(exch::datayx[i]);}
                    CUDA_CALL(cudaMemcpy(Cdyx,otemp,exch::datayx.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datayz.size() ; i++){otemp[i]=static_cast<double>(exch::datayz[i]);}
                    CUDA_CALL(cudaMemcpy(Cdyz,otemp,exch::datayz.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datazx.size() ; i++){otemp[i]=static_cast<double>(exch::datazx[i]);}
                    CUDA_CALL(cudaMemcpy(Cdzx,otemp,exch::datazx.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datazy.size() ; i++){otemp[i]=static_cast<double>(exch::datazy[i]);}
                    CUDA_CALL(cudaMemcpy(Cdzy,otemp,exch::datazy.size()*sizeof(double),cudaMemcpyHostToDevice));
                }
            }
            else if(config::exchm==2)//CSR
            {
                CUDA_CALL(cudaMalloc((void**)&Cxadj,exch::xadj.size()*sizeof(unsigned int)));
                CUDA_CALL(cudaMalloc((void**)&Cadjncy,exch::adjncy.size()*sizeof(unsigned int)));
                CUDA_CALL(cudaMemcpy(Cxadj,exch::xadj.ptr(),exch::xadj.size()*sizeof(unsigned int),cudaMemcpyHostToDevice));
                CUDA_CALL(cudaMemcpy(Cadjncy,exch::adjncy.ptr(),exch::adjncy.size()*sizeof(unsigned int),cudaMemcpyHostToDevice));
                CUDA_CALL(cudaMalloc((void**)&Cdxx,exch::dataxx.size()*sizeof(double)));
                CUDA_CALL(cudaMalloc((void**)&Cdyy,exch::datayy.size()*sizeof(double)));
                CUDA_CALL(cudaMalloc((void**)&Cdzz,exch::datazz.size()*sizeof(double)));
                //dataxx, datayy and datazz should all be the same size
                double *temp=new double[exch::dataxx.size()];
                for(unsigned int i = 0 ; i < exch::dataxx.size() ; i++)
                {
                    temp[i]=static_cast<double>(exch::dataxx[i]);
                }
                CUDA_CALL(cudaMemcpy(Cdxx,temp,exch::dataxx.size()*sizeof(double),cudaMemcpyHostToDevice));
                if(exch::dataxx.size()!=exch::datayy.size())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Something has gone wrong, the length of the diagonal (xx,yy and zz) components of the exchange\ntensor in the CSR format (exch::dataxx etc) should be the same. Here xx and yy not equal.");
                }
                for(unsigned int i = 0 ; i < exch::dataxx.size() ; i++)
                {
                    temp[i]=static_cast<double>(exch::datayy[i]);
                }
                CUDA_CALL(cudaMemcpy(Cdyy,temp,exch::datayy.size()*sizeof(double),cudaMemcpyHostToDevice));
                if(exch::dataxx.size()!=exch::datazz.size())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Something has gone wrong, the length of the diagonal (xx,yy and zz) components of the exchange\ntensor in the CSR format (exch::dataxx etc) should be the same. Here xx and yy not equal.");
                }
                for(unsigned int i = 0 ; i < exch::dataxx.size() ; i++)
                {
                    temp[i]=static_cast<double>(exch::datayy[i]);
                }
                CUDA_CALL(cudaMemcpy(Cdzz,temp,exch::datazz.size()*sizeof(double),cudaMemcpyHostToDevice));
                if(config::offdiag==true)
                {
                    unsigned int arsize=exch::dataxy.size();
                    if(!(exch::dataxy.size()==arsize && exch::dataxz.size()==arsize && exch::datayx.size()==arsize && exch::datayz.size()==arsize && exch::datazx.size()==arsize && exch::datazy.size()==arsize))
                    {
                        error::errPreamble(__FILE__,__LINE__);
                    //    std::cout << exch::dataxy.size() << "\t"<< exch::dataxz.size() << "\t"<< exch::datayx.size() << "\t"<< exch::datayz.size() << "\t"<< exch::datazx.size() << "\t"<< exch::datazy.size() << std::endl;
                        error::errMessage("Size of the off diagonal exchange components (in DIA format) are not of the same length.");
                    }
                    CUDA_CALL(cudaMalloc((void**)&Cdxy,exch::dataxy.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdxz,exch::dataxz.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdyx,exch::datayx.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdyz,exch::datayz.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdzx,exch::datazx.size()*sizeof(double)));
                    CUDA_CALL(cudaMalloc((void**)&Cdzy,exch::datazy.size()*sizeof(double)));
                    double *otemp=new double[exch::dataxy.size()];for(unsigned int i = 0 ; i < exch::dataxy.size() ; i++){otemp[i]=static_cast<double>(exch::dataxy[i]);}
                    CUDA_CALL(cudaMemcpy(Cdxy,otemp,exch::dataxy.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::dataxz.size() ; i++){otemp[i]=static_cast<double>(exch::dataxz[i]);}
                    CUDA_CALL(cudaMemcpy(Cdxz,otemp,exch::dataxz.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datayx.size() ; i++){otemp[i]=static_cast<double>(exch::datayx[i]);}
                    CUDA_CALL(cudaMemcpy(Cdyx,otemp,exch::datayx.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datayz.size() ; i++){otemp[i]=static_cast<double>(exch::datayz[i]);}
                    CUDA_CALL(cudaMemcpy(Cdyz,otemp,exch::datayz.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datazx.size() ; i++){otemp[i]=static_cast<double>(exch::datazx[i]);}
                    CUDA_CALL(cudaMemcpy(Cdzx,otemp,exch::datazx.size()*sizeof(double),cudaMemcpyHostToDevice));
                    for(unsigned int i = 0 ; i < exch::datazy.size() ; i++){otemp[i]=static_cast<double>(exch::datazy[i]);}
                    CUDA_CALL(cudaMemcpy(Cdzy,otemp,exch::datazy.size()*sizeof(double),cudaMemcpyHostToDevice));
                }
            }
        }

        if(exch::inc4spin)
        {
            CUDA_CALL(cudaMalloc((void**)&Cxadj_jkl,exch::xadj_j.size()*sizeof(unsigned int)));
            CUDA_CALL(cudaMalloc((void**)&Cadjncy_j,exch::adjncy_j.size()*sizeof(unsigned int)));
            CUDA_CALL(cudaMalloc((void**)&Cadjncy_k,exch::adjncy_k.size()*sizeof(unsigned int)));
            CUDA_CALL(cudaMalloc((void**)&Cadjncy_l,exch::adjncy_l.size()*sizeof(unsigned int)));
            CUDA_CALL(cudaMemcpy(Cxadj_jkl,exch::xadj_j.ptr(),exch::xadj_j.size()*sizeof(unsigned int),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(Cadjncy_j,exch::adjncy_j.ptr(),exch::adjncy_j.size()*sizeof(unsigned int),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(Cadjncy_k,exch::adjncy_k.ptr(),exch::adjncy_k.size()*sizeof(unsigned int),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(Cadjncy_l,exch::adjncy_l.ptr(),exch::adjncy_l.size()*sizeof(unsigned int),cudaMemcpyHostToDevice));
        }
        CUDA_CALL(cudaMalloc((void**)&CHDemag,3*geom::nspins*sizeof(float)));
        CUDA_CALL(cudaMalloc((void**)&CHstg,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&CInitHstg,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Cspin,3*geom::nspins*sizeof(double)));
        if(llg::intscheme==0)
        {
            CUDA_CALL(cudaMalloc((void**)&Cespin,3*geom::nspins*sizeof(double)));
            CUDA_CALL(cudaMalloc((void**)&Cfn,3*geom::nspins*sizeof(double)));
        }
        else if(llg::intscheme==1)//RK4
        {
            CUDA_CALL(cudaMalloc((void**)&CRK4k1,3*geom::nspins*sizeof(double)));
            CUDA_CALL(cudaMalloc((void**)&CRK4k2,3*geom::nspins*sizeof(double)));
            CUDA_CALL(cudaMalloc((void**)&CRK4k3,3*geom::nspins*sizeof(double)));
            CUDA_CALL(cudaMalloc((void**)&CRK4k4,3*geom::nspins*sizeof(double)));
        }

        CUDA_CALL(cudaMalloc((void**)&CDetFields,3*geom::nspins*sizeof(double)));
        curandN=3*geom::nspins;
        if(curandN%2!=0)//Check if we have an odd number of random numbers (MUST BE EVEN!!!)
        {
            curandN=curandN+1;
        }
        CUDA_CALL(cudaMalloc((void**)&Crand,curandN*sizeof(float)));
        CUDA_CALL(cudaMalloc((void**)&Ck1udir,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&CH,3*geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Clambda,geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Csigma,geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Cllgpf,geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Ck1u,geom::nspins*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Ckx,geom::nspins*sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc((void**)&Cky,geom::nspins*sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc((void**)&Ckz,geom::nspins*sizeof(unsigned int)));
        //--------------------------------------------------------------------------------
        //this section sorts out the copying of the data from the CPU to the card
        //--------------------------------------------------------------------------------
        //copy the sigma prefactor
        CUDA_CALL(cudaMemcpy(Csigma,geom::sigma.ptr(),geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Clambda,geom::lambda.ptr(),geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Cllgpf,geom::llgpf.ptr(),geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Ck1u,anis::k1u.ptr(),geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
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
        //copy the 3 Nspin long arrays containing the staggered field to the tnsfa array
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            tnsda[3*i]=  fields::Hstagx(i);
            tnsda[3*i+1]=fields::Hstagy(i);
            tnsda[3*i+2]=fields::Hstagz(i);
        }
        CUDA_CALL(cudaMemcpy(CHstg,tnsda,geom::nspins*3*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(CInitHstg,tnsda,geom::nspins*3*sizeof(double),cudaMemcpyHostToDevice));

        //copy the first order uniaxial anisotropy direction
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            for(unsigned int coord = 0 ; coord < 3 ; coord++)
            {
                tnsda[3*i+coord]=anis::k1udir(i,coord);
            }
        }
        CUDA_CALL(cudaMemcpy(Ck1udir,tnsda,geom::nspins*3*sizeof(double),cudaMemcpyHostToDevice));
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
        util::copy3vecto1(geom::nspins,spins::Sx,spins::Sy,spins::Sz,tnsda);
        //copy spin data to single array
        //copy spin data to card
        FIXOUT(config::Log,"Copying spin memory to device:" << std::flush);
        CUDA_CALL(cudaMemcpy(Cspin,tnsda,3*geom::nspins*sizeof(double),cudaMemcpyHostToDevice));
        SUCCESS(config::Log);
        //zero the field array
        for(unsigned int i = 0 ; i < 3*geom::nspins ; i++){tnsfa[i]=0.0;}CUDA_CALL(cudaMemcpy(CH,tnsfa,3*geom::nspins*sizeof(float),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(CHDemag,tnsfa,3*geom::nspins*sizeof(float),cudaMemcpyHostToDevice));

        if(config::exchm==0)
        {
            //and copy the species list
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                nsia[i]=geom::lu(i,3);
            }
            CUDA_CALL(cudaMemcpy(Cspec,nsia,geom::nspins*sizeof(unsigned int),cudaMemcpyHostToDevice));
        }
        else if(config::exchm>0 && config::dipm==0 && config::inc_dip==true)
        {
            //and copy the species list
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                nsfa[i]=geom::mu[i];
            }
            CUDA_CALL(cudaMemcpy(Cmagmom,nsfa,geom::nspins*sizeof(float),cudaMemcpyHostToDevice));
        }


        //make sure we clean up when the program exits
        atexit(deallocate_cuda_memory);
    }
    void CsetStagFields()
    {
        double *tnsda=new double[3*geom::nspins];
        //copy the 3 Nspin long arrays containing the staggered field to the tnsfa array
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            tnsda[3*i]=  fields::Hstagx(i);
            tnsda[3*i+1]=fields::Hstagy(i);
            tnsda[3*i+2]=fields::Hstagz(i);
        }
        CUDA_CALL(cudaMemcpy(CHstg,tnsda,geom::nspins*3*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(CInitHstg,tnsda,geom::nspins*3*sizeof(double),cudaMemcpyHostToDevice));
    }
    void CsetStagFieldsZero()
    {
        double *tnsda=new double[3*geom::nspins];
        //copy the 3 Nspin long arrays containing the staggered field to the tnsfa array
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            tnsda[3*i]=  0.0;
            tnsda[3*i+1]=0.0;
            tnsda[3*i+2]=0.0;
        }
        CUDA_CALL(cudaMemcpy(CHstg,tnsda,geom::nspins*3*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(CInitHstg,tnsda,geom::nspins*3*sizeof(double),cudaMemcpyHostToDevice));
    }
}
