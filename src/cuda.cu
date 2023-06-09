// File: cuda.cu
// Author:Tom Ostler
// Last-modified: 19 May 2023 13:56:33
// Formerly cuLLB.cu
#include "../inc/cuda.h"
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
#include "../inc/cuintRK4.h"
#include "../inc/cuintRK5.h"
#include "../inc/llg.h"
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
namespace cullg
{


    void rampStag(double rampval)
    {
        cufields::CincStagFields<<<blockspergrid,threadsperblock>>>(geom::nspins,rampval,CInitHstg,CHstg);
    }

    void llgGPU(unsigned int& t)
    {
        //in this case we are using the interaction matrix for both the exchange and the
        //dipole-dipole field so we might aswell update both at once
        if(config::exchm==0)
        {
            cufields::CZero5DRSArrays<<<rsarzpblockspergrid,threadsperblock>>>(geom::zps*3*geom::ucm.GetNMS(),CHr,CSr,CHk,CSk);
            //copy the spin data to the zero padded arrays
            cufields::CCopySpin<<<blockspergrid,threadsperblock>>>(geom::nspins,Cspin,CSr,Ckx,Cky,Ckz,Cspec);
            //forward transform
            spins_forward();
            //perform convolution
            cufields::CFConv<<<zpblockspergrid,threadsperblock>>>(geom::zps,geom::ucm.GetNMS(),CNk,CHk,CSk);
            //transform the fields back
            fields_back();
            //copy the fields from the zero padded array to the demag field array
            cufields::CCopyFields<<<blockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CH,CHr,Ckx,Cky,Ckz,Cspec);
        }
        else if(config::exchm>0 && config::dipm==0 && config::inc_dip==true)
        {
            //calculate the demag field once per spins::update steps and store in CHDemag. This way when we do the
            //matrix multiplication we just add the exchange field to the demag part.
            if(t%spins::update==0)
            {
                cufields::CZero4DRSArrays<<<rsarzpblockspergrid,threadsperblock>>>(geom::zps*3,CHr,CSr,CHk,CSk);
                //copy the spin data to the zero padded arrays
                cufields::CdipCopySpin<<<blockspergrid,threadsperblock>>>(geom::nspins,Cspin,CSr,Ckx,Cky,Ckz,Cmagmom);
                //copy the fields from the zero padded array to the demag field array
//                cufields::CdipCopyFields<<<zpblockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CHDemag,CHr,Ckx,Cky,Ckz);
                //forward transform
                spins_forward();
                cufields::CdipFConv<<<zpblockspergrid,threadsperblock>>>(geom::zps,CNk,CHk,CSk);
                fields_back();
                cufields::CdipCopyFields<<<blockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CHDemag,CHr,Ckx,Cky,Ckz);
                
                //copy spin arrays back to CPU
                //FOR DEBUGGING THE DIPOLAR FIELD
                    /*float *temp=NULL;
                    temp = new float [3*geom::nspins];
                    CUDA_CALL(cudaMemcpy(temp,CHDemag,3*geom::nspins*sizeof(float),cudaMemcpyDeviceToHost));
                    for(unsigned int i = 0 ; i < geom::nspins ; i++)
                    {
                        std::cout << i << "\t" << temp[3*i] << "\t" << temp[3*i+1] << "\t" << temp[3*i+2] << std::endl;
                    }
                    exit(0);*/
                

            }
            
            if(config::exchm==1)//DIA
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cspin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                }
            }
            else if(config::exchm==2)//CSR
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cspin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                }
            }

        }
        else if(config::exchm>0 && config::inc_dip==false)
        {

            if(config::exchm==1)//DIA
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cspin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                }
            }
            else if(config::exchm==2)//CSR
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cspin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                }
            }
        }
        //calcalute the four spin term?
        if(exch::inc4spin)
        {
            cufields::CSpMV_CSR_FourSpin<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj_jkl,Cadjncy_j,Cadjncy_k,Cadjncy_l,CH,Cspin);//,CH,Cspin);
        }
        //generate the random numbers
        CURAND_CALL(curandGenerateNormal(gen,Crand,curandN,0.0,1.0));
/*            float *temp=NULL;
            temp = new float [3*geom::nspins];
            CUDA_CALL(cudaMemcpy(temp,CH,3*geom::nspins*sizeof(float),cudaMemcpyDeviceToHost));
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                std::cout << std::setprecision(10) << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << temp[3*i] << "\t" << temp[3*i+1] << "\t" << temp[3*i+2] << std::endl;
            }
            delete [] temp;
            temp=NULL;
            exit(0);*/
        cuint::CHeun1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,Crand,Cfn,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
        //in this case we are using the interaction matrix for both the exchange and the
        //dipole-dipole field so we might aswell update both at once
        if(config::exchm==0)
        {
            cufields::CZero5DRSArrays<<<rsarzpblockspergrid,threadsperblock>>>(geom::zps*3*geom::ucm.GetNMS(),CHr,CSr,CHk,CSk);
            //copy the spin data to the zero padded arrays
            cufields::CCopySpin<<<blockspergrid,threadsperblock>>>(geom::nspins,Cespin,CSr,Ckx,Cky,Ckz,Cspec);
            //forward transform
            spins_forward();
            //perform convolution
            cufields::CFConv<<<zpblockspergrid,threadsperblock>>>(geom::zps,geom::ucm.GetNMS(),CNk,CHk,CSk);
            //transform the fields back
            fields_back();
            //copy the fields from the zero padded array to the demag field array
            cufields::CCopyFields<<<blockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CH,CHr,Ckx,Cky,Ckz,Cspec);
        }
        else if(config::exchm>0 && config::dipm==0 && config::inc_dip==true)
        {
            //we don't want to update the dipole field here
            if(config::exchm==1)//DIA
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                }
            }
            else if(config::exchm==2)//CSR
            {
                if(config::offdiag)
                {
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                }
            }

        }
        else if(config::exchm>0 && config::inc_dip==false)
        {
            if(config::exchm==1)//DIA
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                }
            }
            else if(config::exchm==2)//CSR
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                }
            }
        }
        //calcalute the four spin term?
        if(exch::inc4spin)
        {
            cufields::CSpMV_CSR_FourSpin<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj_jkl,Cadjncy_j,Cadjncy_k,Cadjncy_l,CH,Cespin);
        }
        cuint::CHeun2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,Crand,Cfn,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CDetFields,CHstg);
        if(t%spins::update==0)
        {
            //copy spin arrays back to CPU
            double *temp=NULL;
            temp = new double [3*geom::nspins];
            CUDA_CALL(cudaMemcpy(temp,Cspin,3*geom::nspins*sizeof(double),cudaMemcpyDeviceToHost));
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                spins::Sx[i]=temp[3*i];
                spins::Sy[i]=temp[3*i+1];
                spins::Sz[i]=temp[3*i+2];
//                				std::cout << spins::Sx[i] << "\t" << spins::Sy[i] << "\t" << spins::Sz[i] << "\t" << sqrt(spins::Sx[i]*spins::Sx[i] + spins::Sy[i]*spins::Sy[i] + spins::Sz[i]*spins::Sz[i])<< std::endl;
            }
            CUDA_CALL(cudaMemcpy(temp,CDetFields,3*geom::nspins*sizeof(double),cudaMemcpyDeviceToHost));
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                fields::Hx[i]=temp[3*i];
                fields::Hy[i]=temp[3*i+1];
                fields::Hz[i]=temp[3*i+2];
            }
            delete [] temp;
            temp=NULL;
        }

    }


    void cuinit()
    {

        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**CUDA details***" << std::endl;
        FIXOUT(config::Info,"Resetting device:" << std::flush);
        CUDA_CALL(cudaDeviceReset());
        SUCCESS(config::Info);

        if(!config::cfg.exists("cuda"))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting \"cuda\" does not exist. Check your config file.");
        }
        libconfig::Setting &setting = config::cfg.lookup("cuda");
        config::Info << std::noshowpos;
        //FIXOUT(config::Info,"NVCC Compiler:" << COMP << std::endl);
        int device_count=0;
        int device=0;
        //---------------------------------------------------------------
        //Get some of the device properties
        if((cudaGetDevice(&device))!=cudaSuccess)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not get device id");
        }
        if((cudaGetDeviceCount(&device_count))!=cudaSuccess)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not get number of devices");
        }
        if(device>device_count)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("GPU device greater than count of devices.");
        }
        if(cudaGetDeviceProperties(&deviceProp,device)!=cudaSuccess)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not get device properties");
        }
        FIXOUT(config::Info,"Selecting device:" << std::endl);
        for(unsigned int i = 0 ; i < device_count ; i++)
        {
            std::stringstream sstr;
            sstr << "Setting device " << i;
            std::string str=sstr.str();
            FIXOUT(config::Info,str << std::flush);
            if(cudaSetDevice(device)!=cudaSuccess)
            {
                config::Info << "Failed" << std::endl;
            }
            else
            {
                config::Info << "Success" << std::endl;
            }
        }
        FIXOUT(config::Info,"Number of devices:" << device_count << std::endl);
        FIXOUT(config::Info,"Device selected:" << device << std::endl);
        FIXOUT(config::Info,"Device major.minor:" << deviceProp.major << "." << deviceProp.minor << std::endl);
        //---------------------------------------------------------------
        setting.lookupValue("threadsperblock",threadsperblock);
        FIXOUT(config::Info,"Number of threads per block:" << threadsperblock << std::endl);
        blockspergrid=(geom::nspins+threadsperblock-1)/threadsperblock;
        zpblockspergrid=(geom::zps+threadsperblock-1)/threadsperblock;
        //This is the number of block per grid for addressing the elements of the real space
        //spin and field arrays (dimensions: NUMSPEC x 3 x ZPDIM[0] x ZPDIM[1] x ZPDIM[2]
        if(config::exchm==0)
        {
            rsarzpblockspergrid=(geom::zps*3*geom::ucm.GetNMS()+threadsperblock-1)/threadsperblock;
        }
        else if(config::exchm>0)
        {
            rsarzpblockspergrid=(geom::zps*3+threadsperblock-1)/threadsperblock;
        }
        //Same as the rsarzpblockspergrid but with the ZPDIM[2] dimension now replaced with the (ZPDIM[2]+1)/2
        FIXOUT(config::Info,"Blocks per grid:" << blockspergrid << std::endl);
        FIXOUT(config::Info,"Blocks per grid for zero pad workspace:" << zpblockspergrid << std::endl);
        FIXOUT(config::Info,"Blocks per grid for addressing each 5D array:" << rsarzpblockspergrid << std::endl);
        FIXOUT(config::Info,"Device maximum threads per block:" << deviceProp.maxThreadsPerBlock << std::endl);
        FIXOUT(config::Info,"Device registers per block:" << deviceProp.regsPerBlock << std::endl);
        FIXOUT(config::Info,"Device total const memory:" << deviceProp.totalConstMem << " (bytes)" << std::endl);
        FIXOUT(config::Info,"Device total global memory:" << deviceProp.totalGlobalMem << " (bytes)" << std::endl);

        unsigned long long int curandseed=config::seed;
        FIXOUT(config::Info,"Curand seed:" << curandseed << std::endl);
        //initialize the random number generator
        check_cuda_errors(__FILE__,__LINE__);
        FIXOUT(config::Info,"Initializing curand random number generator" << std::flush);
        if((curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_DEFAULT))!=CURAND_STATUS_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CURAND failed to create random number generator");
        }
        check_cuda_errors(__FILE__,__LINE__);
        if((curandSetPseudoRandomGeneratorSeed(gen,curandseed))!=CURAND_STATUS_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CURAND failed to set random number seed");
        }
        check_cuda_errors(__FILE__,__LINE__);

        if((curandGenerateSeeds(gen))!=CURAND_STATUS_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CURAND failed to generate random number generator seeds");
        }
        check_cuda_errors(__FILE__,__LINE__);
        /*if((cudaThreadSetLimit(cudaLimitStackSize,1024))!=cudaSuccess)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CUDA ERROR: Failed to set thread limit");
        }*/

        config::Info << "Done" << std::endl;
        FIXOUT(config::Info,"Checking for any cuda errors:" << std::flush);
        check_cuda_errors(__FILE__,__LINE__);
        config::Info << "Done" << std::endl;

        FIXOUT(config::Info,"Allocating memory on device" << std::flush);
        allocate_memory_on_card();
        config::Info << "Done" << std::endl;
        FIXOUT(config::Info,"Copying fourier transformed interaction matrix to device:" << std::flush);
        setup_fourier_transform();
        config::Info << "Done" << std::endl;
        config::printline(config::Info);
        //__constant__ memory only have .cu scope therefore to use the variables 
        //the variables have to be declared in each .cu file and the variables initialized.
        if(llg::intscheme==0)
        {
            cuint::copyConstData();
        }
        else if(llg::intscheme==1)
        {
            cuintRK4::copyConstData();
        }
        else if(llg::intscheme==2)
        {
            cuintRK5::copyConstData();
        }
        cufields::copyConstData();
    }
}
