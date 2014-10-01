// File: cuda.cu
// Author:Tom Ostler
// Last-modified: 01 Oct 2014 18:24:38
// Formerly cuLLB.cu
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

    void initGPU()
    {
        CUDA_CALL(cudaDeviceReset());
    }
    void llgGPU(unsigned int& t)
    {

        //copy the spin data to the zero padded arrays
        cufields::CCopySpin<<<zpblockspergrid,threadsperblock>>>(geom::nspins,Cspin,CSr,Ckx,Cky,Ckz,Cspec);
        //forward transform
        spins_forward();
        //perform convolution
        cufields::CFConv<<<czpblockspergrid,threadsperblock>>>(geom::czps,geom::ucm.GetNMS(),CNk,CHk,CSk);
        //transform the fields back
        fields_back();
        //copy the fields from the zero padded array to the demag field array
        cufields::CCopyFields<<<blockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CH,CHr,Ckx,Cky,Ckz,Cspec);
        //FOR DEBUGGING THE DIPOLAR FIELD/
        float temp1[3*geom::nspins];
        CUDA_CALL(cudaMemcpy(temp1,CH,3*geom::nspins*sizeof(float),cudaMemcpyDeviceToHost));
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            int ijk[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            std::cout << i << "\t" << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << "\t" << temp1[3*i] << "\t" << temp1[3*i+1] << "\t" << temp1[3*i+2] << std::endl;

        }
        exit(0);

        //generate the random numbers
        CURAND_CALL(curandGenerateNormal(gen,Crand,3*geom::nspins,0.0,1.0));
        cuint::CHeun1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,Crand,Cfn,Csigma,Cllgpf,Clambda);
        cufields::CCopySpin<<<zpblockspergrid,threadsperblock>>>(geom::nspins,Cspin,CSr,Ckx,Cky,Ckz,Cspec);
        //forward transform
        spins_forward();
        //perform convolution
        cufields::CFConv<<<czpblockspergrid,threadsperblock>>>(geom::czps,geom::ucm.GetNMS(),CNk,CHk,CSk);
        //transform the fields back
        fields_back();
        //copy the fields from the zero padded array to the demag field array
        cufields::CCopyFields<<<blockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CH,CHr,Ckx,Cky,Ckz,Cspec);

        cuint::CHeun2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,Crand,Cfn,Csigma,Cllgpf,Clambda);
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
                //				std::cout << spins::Sx[i] << "\t" << spins::Sy[i] << "\t" << spins::Sz[i] << "\t" << sqrt(spins::Sx[i]*spins::Sx[i] + spins::Sy[i]*spins::Sy[i] + spins::Sz[i]*spins::Sz[i])<< std::endl;
            }
            delete [] temp;
            temp=NULL;
        }

    }


    void cuinit(int argc,char *argv[])
    {

        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**CUDA details***" << std::endl;
        FIXOUT(config::Info,"Resetting device:" << std::flush);
        CUDA_CALL(cudaDeviceReset());
        SUCCESS(config::Info);

        //the rank of the fourier transform
        nrank=3;
        try
        {
            config::cfg.readFile(argv[1]);
        }
        catch(const libconfig::FileIOException &fioex)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("I/O error while reading config file");
        }
        catch(const libconfig::ParseException &pex)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::cerr << ". Parse error at " << pex.getFile()  << ":" << pex.getLine() << "-" << pex.getError() << "***\n" << std::endl;
            exit(EXIT_FAILURE);
        }
        libconfig::Setting &setting = config::cfg.lookup("cuda");
        config::Info << std::noshowpos;
        FIXOUT(config::Info,"NVCC Compiler:" << COMP << std::endl);
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
        FIXOUT(config::Info,"Number of devices:" << device_count << std::endl);
        FIXOUT(config::Info,"Device selected:" << device << std::endl);
        FIXOUT(config::Info,"Device major.minor:" << deviceProp.major << "." << deviceProp.minor << std::endl);
        //---------------------------------------------------------------
        setting.lookupValue("threadsperblock",threadsperblock);
        FIXOUT(config::Info,"Number of threads per block:" << threadsperblock << std::endl);
        blockspergrid=(geom::nspins+threadsperblock-1)/threadsperblock;
        zpblockspergrid=(geom::zps+threadsperblock-1)/threadsperblock;
        czpblockspergrid=(geom::czps+threadsperblock-1)/threadsperblock;
        FIXOUT(config::Info,"Blocks per grid:" << blockspergrid << std::endl);
        FIXOUT(config::Info,"Blocks per grid for zero pad workspace:" << zpblockspergrid << std::endl);
        FIXOUT(config::Info,"Blocks per grid for complex zero pad workspace:" << czpblockspergrid << std::endl);
        FIXOUT(config::Info,"Device maximum threads per block:" << deviceProp.maxThreadsPerBlock << std::endl);
        FIXOUT(config::Info,"Device registers per block:" << deviceProp.regsPerBlock << std::endl);
        FIXOUT(config::Info,"Device total const memory:" << deviceProp.totalConstMem << " (bytes)" << std::endl);
        FIXOUT(config::Info,"Device total global memory:" << deviceProp.totalGlobalMem << " (bytes)" << std::endl);

        unsigned long long int curandseed=config::seed;
        FIXOUT(config::Info,"Curand seed:" << curandseed << std::endl);
        //initialize the random number generator
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
        config::Info << "NVIDIA-SMI output:\n" << util::exec("nvidia-smi");
        //__constant__ memory only have .cu scope therefore to use the variables 
        //the variables have to be declared in each .cu file and the variables initialized.
        cuint::copyConstData();
        cufields::copyConstData();
    }
}
