// File: cuda.cu
// Author:Tom Ostler
// Last-modified: 10 Jan 2013 13:47:49
// Formally cuLLB.cu
#include "../inc/cufields.h"
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/sim.h"
#include "../inc/intmat.h"
#include "../inc/util.h"
#include "../inc/neigh.h"
#include "../inc/fields.h"
#include "../inc/tdp.h"
#include "../inc/cuint.h"
#include "../inc/array3d.h"
#include "../inc/array.h"
#include "../inc/mat.h"
#include "../inc/mf.h"
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
#ifdef DEBUG
#define CUDA_CALL(x) do { if((x) != cudaSuccess) {\
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CUDA_CALL(x) (x)
#endif /*DEBUG*/
#ifdef DEBUG
#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CURAND_CALL(x) (x)
#endif /*DEBUG*/
#ifdef DEBUG
#define CUFFT_CALL(x) do { if((x) != CUFFT_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CUFFT_CALL(x) (x)
#endif /*DEBUG*/
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace cullb
{
    cudaDeviceProp deviceProp;
    curandGenerator_t gen;
    bool cullbi=false;
    bool cullbmi=false;
    bool cullbffti=false;
    bool cullbinit=false;
    //number of threads per block and blocks per grid
    int threadsperblock,blockspergrid;
    //same but for the zero padded work spaces
    int zpblockspergrid;
    //rank of the FFT
    int nrank=3;
    //device pointers for Fourier space calculations
    static  cufftComplex *CCNxx=NULL;
    static  cufftComplex *CCNxy=NULL;
    static  cufftComplex *CCNxz=NULL;
    static  cufftComplex *CCNyx=NULL;
    static  cufftComplex *CCNyy=NULL;
    static  cufftComplex *CCNyz=NULL;
    static  cufftComplex *CCNzx=NULL;
    static  cufftComplex *CCNzy=NULL;
    static  cufftComplex *CCNzz=NULL;
    static  cufftComplex *CCSx=NULL;
    static  cufftComplex *CCSy=NULL;
    static  cufftComplex *CCSz=NULL;
    static  cufftComplex *CCHx=NULL;
    static  cufftComplex *CCHy=NULL;
    static  cufftComplex *CCHz=NULL;
    //device pointers
    static  double *Cspin=NULL;
    static  double *Cespin=NULL;
    static  int *Cxadj=NULL;
    static  int *Cadjncy=NULL;
    static  float *Cfspin=NULL;
    static  float *CTemp=NULL;
    static  float *Crand=NULL;
    static  float *Ccoord=NULL;
    static  float *CHDemag=NULL;
    static  int *Clu=NULL;
    static  int *Czpsn=NULL;//The is the zero pad spin number
	static  double *CsurfArea=NULL;
	static  double *CW1pf=NULL;
	static  double *CW2pf=NULL;
	static  double *Csigma=NULL;
    static  double *Cfn=NULL;
    //cufft plans
    cufftHandle C3DP;

    void initGPU()
    {
        //CUDA_CALL(cudaDeviceReset());
        float *nsfa=new float[geom::ss];
        util::fillfloat(geom::ss,tdp::systemp.ptr(),nsfa);
        CUDA_CALL(cudaMemcpy(CTemp,nsfa,geom::ss*sizeof(float),cudaMemcpyHostToDevice));
        double dtemp[3*geom::ss];
        util::copy3vecto1(geom::ss,spins::sx,spins::sy,spins::sz,dtemp);
        CUDA_CALL(cudaMemcpy(Cspin,dtemp,3*geom::ss*sizeof(double),cudaMemcpyHostToDevice));
        //copy the contents of the Mean Field arrays to the device
        cuint::copyConstData();
        //allocate the function pointers used to calculate temperature
        //dependent properties on the card
        int lmetype=0,lchipartype=0,lchiperptype=0,lexchstifftype=0;
        int lmefunc=-1,lchiparfunc=-1,lchiperpfunc=-1,lexchstifffunc=-1;
        if(tdp::metype=="function")
        {
            lmetype=0;
            lmefunc=tdp::mefunc;
            if(lmefunc!=0)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Magnetization function type not currently supported on GPU");
            }
        }
        else if(tdp::metype=="meanfield")
        {
            lmetype=-1;
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Magnetization calculation not supported on GPU");
        }
        if(tdp::chipartype=="function")
        {
            lchipartype=0;
            lchiparfunc=tdp::chiparfunc;
            if(lchiparfunc!=0)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Chi par function type not currently supported on GPU");
            }

        }
        else if(tdp::chipartype=="meanfield")
        {
            lchipartype=-1;

        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Parallel susceptibility calculation not supported on GPU");
        }
        if(tdp::chiperptype=="function")
        {
            lchiperptype=0;
            lchiperpfunc=tdp::chiperpfunc;
            if(lchiperpfunc!=0)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Chi perp function type not currently supported on GPU");
            }
        }
        else if(tdp::chiperptype=="meanfield")
        {
            lchiperptype=-1;

        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Perpendicular susceptibility calculation not supported on GPU");
        }
		if(fields::checkexchange==true && tdp::exchstifftype=="function")
		{
			lexchstifftype=0;
			lexchstifffunc=tdp::exchstifffunc;
			if(lexchstifffunc!=0)
			{
				error::errPreamble(__FILE__,__LINE__);
				error::errMessage("Exchange stiffness function type not currently supported on GPU");
			}
		}
        else if(fields::checkexchange==false)
        {
            lexchstifftype=-1;
        }
		else
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Exchange stiffness calculation type not supported on GPU");
		}

        cuint::CSetFunctionPointer<<<threadsperblock,blockspergrid>>>(lmetype,lmefunc,lchipartype,lchiparfunc,lchiperptype,lchiperpfunc,lexchstifftype,lexchstifffunc,fields::anist,fields::ff);

        cullbinit=true;

    }
    void LLBGPU(unsigned int& t)
    {

        assert(geom::gi);
        assert(config::lcf);
        assert(neigh::ni);
        assert(util::ui);
        assert(cullbi);
        assert(cullbmi);
        assert(fields::fi);
        assert(tdp::tdpi);
        if(fields::checkdipolar)
        {
            if(t==0 || t%fields::dfu==0)
            {
                if(fields::dfc=="fft")
                {

                    assert(cullbffti);
                    //copy the spin data to the zero padded arrays
                    cufields::CCopySpin<<<threadsperblock,zpblockspergrid>>>(geom::zps,geom::ss,Cspin,Czpsn,CCSx,CCSy,CCSz);
                    //forward transform
                    spins_forward();
                    //perform convolution
                    cufields::CFConv<<<threadsperblock,zpblockspergrid>>>(geom::zps,CCNxx,CCNxy,CCNxz,CCNyx,CCNyy,CCNyz,CCNzx,CCNzy,CCNzz,CCHx,CCHy,CCHz,CCSx,CCSy,CCSz);
                    //transform the fields back
                    fields_back();
                    //copy the fields from the zero padded array to the demag field array
                    cufields::CCopyFields<<<threadsperblock,blockspergrid>>>(geom::ss,geom::zps,CHDemag,Czpsn,CCHx,CCHy,CCHz);
                }
                else if(fields::dfc=="bruteforce")
                {
                    cufields::CBfDip<<<threadsperblock,blockspergrid>>>(geom::ss,geom::gsV*mat::Ms,CHDemag,Cfspin,Ccoord);
                }
            }
        }
        for(unsigned int i = 0 ; i < fields::dfu ; i++)
        {
            //FOR DEBUGGING THE DIPOLAR FIELD/
            
            float temp[3*geom::ss];
            CUDA_CALL(cudaMemcpy(temp,CHDemag,3*geom::ss*sizeof(float),cudaMemcpyDeviceToHost));
            for(unsigned int i = 0 ; i < geom::ss ; i++)
            {
                int ijk[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                std::cout << i << "\t" << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << "\t" << temp[3*i] << "\t" << temp[3*i+1] << "\t" << temp[3*i+2] << std::endl;
            }
            exit(0);
            
            if(fields::checkthermal==true)
            {
                //generate the random numbers
                CURAND_CALL(curandGenerateNormal(gen,Crand,6*geom::ss,0.0,1.0));
            }
            else
            {
                //or set them to zero. This doesn't need to be called each time.
                cufields::CZeroCrand<<<threadsperblock,blockspergrid>>>(geom::ss,Crand);
            }
            cuint::CHeun1<<<threadsperblock,blockspergrid>>>(geom::ss,float(fields::fH[0]),float(fields::fH[1]),float(fields::fH[2]),mat::Tc,CHDemag,Cfspin,Cspin,Cespin,CTemp,Cxadj,Cadjncy,CsurfArea,Csigma,Crand,Cfn);
/*            check_cuda_errors(__FILE__,__LINE__);
            float temp[geom::ss];
            CUDA_CALL(cudaMemcpy(temp,CTemp,geom::ss*sizeof(float),cudaMemcpyDeviceToHost));
            for(unsigned int i = 0 ; i < geom::ss ; i++)
            {
                std::cout << tdp::systemp[0]+i*10 << "\t" << temp[i] << std::endl;
            }
            exit(0);*/
            cuint::CHeun2<<<threadsperblock,blockspergrid>>>(geom::ss,float(fields::fH[0]),float(fields::fH[1]),float(fields::fH[2]),mat::Tc,Cfspin,Cspin,CHDemag,CTemp,Crand,Cespin,Cfn,Cxadj,Cadjncy,CsurfArea,Csigma);
        }
        t+=fields::dfu;
        double *temp=NULL;
        temp = new double [3*geom::ss];
        CUDA_CALL(cudaMemcpy(temp,Cspin,3*geom::ss*sizeof(double),cudaMemcpyDeviceToHost));
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            spins::sx[i]=temp[3*i];
            spins::sy[i]=temp[3*i+1];
            spins::sz[i]=temp[3*i+2];
        }

    }


    void cuinit(int argc,char *argv[])
    {
        assert(spins::si);
        assert(mat::mi);
        assert(geom::gi);
        assert(config::lcf);
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**CUDA details***" << std::endl;

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

        FIXOUT(config::Info,"NVCC Compiler:" << COMP << std::endl);
        int device_count=0;
        int device=0;
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
        if(deviceProp.major < 2)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Cuda compute capability of 2.0 or greater is required for the use of function pointers.");
        }
        setting.lookupValue("threadsperblock",threadsperblock);
        FIXOUT(config::Info,"Number of threads per block:" << threadsperblock << std::endl);
        blockspergrid=(geom::ss+threadsperblock-1)/threadsperblock;
        zpblockspergrid=(geom::zps+threadsperblock-1)/threadsperblock;
        FIXOUT(config::Info,"Blocks per grid:" << blockspergrid << std::endl);
        FIXOUT(config::Info,"Blocks per grid for zero pad workspace:" << zpblockspergrid << std::endl);
        FIXOUT(config::Info,"Device maximum threads per block:" << deviceProp.maxThreadsPerBlock << std::endl);
        FIXOUT(config::Info,"Device registers per block:" << deviceProp.regsPerBlock << std::endl);
        FIXOUT(config::Info,"Device total const memory:" << deviceProp.totalConstMem << " (bytes)" << std::endl);
        FIXOUT(config::Info,"Device total global memory:" << deviceProp.totalGlobalMem << " (bytes)" << std::endl);
        if((cudaSetDevice(device))!=cudaSuccess)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("cudaSetDevice returned cudaErrorInvalidDevice");
        }
        else
        {
            //double check the device was properly selected
            if(cudaGetDevice(&device)!=cudaSuccess)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not get device on double check");
            }
        }
        //initialize the random number generator
        FIXOUT(config::Info,"Initializing curand random number generator" << std::flush);
        unsigned long long int curandseed=config::seed;
        if((curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_DEFAULT))!=CURAND_STATUS_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CURAND failed to create random number generator");
        }
        if((curandSetPseudoRandomGeneratorSeed(gen,curandseed))!=CURAND_STATUS_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CURAND failed to set random number seed");
        }
        if((curandGenerateSeeds(gen))!=CURAND_STATUS_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CURAND failed to generate random number generator seeds");
        }
        config::Info << "Done" << std::endl;
        FIXOUT(config::Info,"Checking for any cuda errors:" << std::flush);
        check_cuda_errors(__FILE__,__LINE__);
        config::Info << "Done" << std::endl;

        FIXOUT(config::Info,"Allocating memory on device" << std::flush);
        unsigned int memcount=0;
        allocate_memory_on_card(memcount);
        config::Info << "Done" << std::endl;
        FIXOUT(config::Info,"Approximate amount of global memory allocated:" << memcount << " (bytes)" << std::endl);
        if(fields::dfc=="fft")
        {
            FIXOUT(config::Info,"Setting up Fourier transforms on device:" << std::flush);
            setup_fourier_transform();
            config::Info << "Done" << std::endl;
            FIXOUT(config::Info,"Estimated memory for cufft transform:" << double(geom::zps)*double(sizeof(float))*2.0/1024.0/1024.0 << std::endl);
        }
        //size_t f,t;
        //CUcontext context;
        //cuCtxCreate(&context,0,device);
        //CUresult result;

        //result=cuMemGetInfo(&f,&t);
        //FIXOUT(config::Info,"Bytes free on device:" << f << std::endl);
        //FIXOUT(config::Info,"Bytes used on device:" << t << std::endl);
        //This should be asserted before running any kernels
        //so that memory is properly allocated.
        cullbi=true;
    }

    void spins_forward()
    {
        check_cuda_errors(__FILE__,__LINE__);
        assert(cullbmi);
        assert(cullbffti);
        CUFFT_CALL(cufftExecC2C(C3DP,CCSx,CCSx,CUFFT_FORWARD));
        CUFFT_CALL(cufftExecC2C(C3DP,CCSy,CCSy,CUFFT_FORWARD));
        CUFFT_CALL(cufftExecC2C(C3DP,CCSz,CCSz,CUFFT_FORWARD));
    }

    void fields_back()
    {
        assert(cullbmi);
        assert(cullbffti);
        CUFFT_CALL(cufftExecC2C(C3DP,CCHx,CCHx,CUFFT_INVERSE));
        CUFFT_CALL(cufftExecC2C(C3DP,CCHy,CCHy,CUFFT_INVERSE));
        CUFFT_CALL(cufftExecC2C(C3DP,CCHz,CCHz,CUFFT_INVERSE));
    }

    void allocate_memory_on_card(unsigned int& memcount)
    {
        //all of the GPU memory allocations should happen here.
        //--------------------------------------------------------------------------------
        CUDA_CALL(cudaMalloc((void**)&CCNxx,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNxy,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNxz,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNyx,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNyy,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNyz,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNzx,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNzy,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCNzz,geom::zps*sizeof(cufftComplex)));

        CUDA_CALL(cudaMalloc((void**)&CCSx,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCSy,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCSz,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCHx,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCHy,geom::zps*sizeof(cufftComplex)));
        CUDA_CALL(cudaMalloc((void**)&CCHz,geom::zps*sizeof(cufftComplex)));
        memcount+=sizeof(cufftComplex)*geom::zps*12;
        CUDA_CALL(cudaMalloc((void**)&Cspin,3*geom::ss*sizeof(double)));
        CUDA_CALL(cudaMalloc((void**)&Cespin,3*geom::ss*sizeof(double)));
        memcount+=sizeof(double)*3*geom::ss*2;
        CUDA_CALL(cudaMalloc((void**)&Cxadj,neigh::xadj.size()*sizeof(int)));
        memcount+=neigh::xadj.size()*sizeof(int);
        CUDA_CALL(cudaMalloc((void**)&Cadjncy,neigh::adjncy.size()*sizeof(int)));
        memcount+=neigh::adjncy.size()*sizeof(int);
        CUDA_CALL(cudaMalloc((void**)&Cfspin,3*geom::ss*sizeof(float)));
        memcount+=3*geom::ss*sizeof(float);
        CUDA_CALL(cudaMalloc((void**)&Crand,6*geom::ss*sizeof(float)));
        memcount+=6*geom::ss*sizeof(float);
        CUDA_CALL(cudaMalloc((void**)&CTemp,geom::ss*sizeof(float)));
        memcount+=geom::ss*sizeof(float);
        CUDA_CALL(cudaMalloc((void**)&Ccoord,3*geom::ss*sizeof(float)));
        CUDA_CALL(cudaMalloc((void**)&CHDemag,3*geom::ss*sizeof(float)));
        memcount+=3*geom::ss*sizeof(float)*2;
        CUDA_CALL(cudaMalloc((void**)&Clu,3*geom::ss*sizeof(int)));
        memcount+=3*geom::ss*sizeof(int);
        CUDA_CALL(cudaMalloc((void**)&Czpsn,geom::ss*sizeof(int)));
        memcount+=geom::ss*sizeof(int);
		CUDA_CALL(cudaMalloc((void**)&CsurfArea,neigh::surfArea.size()*sizeof(double)));
        memcount+=neigh::surfArea.size()*sizeof(double);
		CUDA_CALL(cudaMalloc((void**)&CW1pf,geom::ss*sizeof(double)));
		CUDA_CALL(cudaMalloc((void**)&CW2pf,geom::ss*sizeof(double)));
		CUDA_CALL(cudaMalloc((void**)&Csigma,geom::ss*sizeof(double)));
        memcount+=geom::ss*sizeof(double);
        CUDA_CALL(cudaMalloc((void**)&Cfn,3*geom::ss*sizeof(double)));
        memcount+=geom::ss*sizeof(double)*3;

        //--------------------------------------------------------------------------------
        //this section sorts out the copying of the data from the CPU to the card
        //--------------------------------------------------------------------------------
        //declare some arrays for doing copying to card
        //Nspins float array, 3*Nspins float array.
        float *nsfa=new float[geom::ss];
        float *tnsfa=new float[3*geom::ss];
        //Nspins double array, 3*Nspins double array
        double *nsda=new double[geom::ss];
        double *tnsda=new double[3*geom::ss];
        //Nspins int array, 3*Nspins int array
        int *nsia=new int[geom::ss];
        int *tnsia=new int[3*geom::ss];
        //copy spin data to single array
        util::copy3vecto1(geom::ss,spins::sx,spins::sy,spins::sz,tnsda);
        //copy spin data to a floating point array
        util::fillfloat(3*geom::ss,tnsda,tnsfa);
        //copy spin data to card
        CUDA_CALL(cudaMemcpy(Cspin,tnsda,3*geom::ss*sizeof(double),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Cfspin,tnsfa,3*geom::ss*sizeof(float),cudaMemcpyHostToDevice));

        //copy temperature data
        util::fillfloat(geom::ss,tdp::systemp.ptr(),nsfa);
        CUDA_CALL(cudaMemcpy(CTemp,nsfa,geom::ss*sizeof(float),cudaMemcpyHostToDevice));
        //copy neighbourlist data
        CUDA_CALL(cudaMemcpy(Cxadj,neigh::xadj.ptr(),neigh::xadj.size()*sizeof(int),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Cadjncy,neigh::adjncy.ptr(),neigh::adjncy.size()*sizeof(int),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CsurfArea,neigh::surfArea.ptr(),neigh::surfArea.size()*sizeof(double),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CW1pf,tdp::sysW1pf.ptr(),tdp::sysW1pf.size()*sizeof(double),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CW2pf,tdp::sysW2pf.ptr(),tdp::sysW2pf.size()*sizeof(double),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(Csigma,tdp::sysSigma.ptr(),geom::ss*sizeof(double),cudaMemcpyHostToDevice));
        
        //take the integer coordinates and calculate the real space coordinates
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                tnsfa[3*i+j]=geom::lu(i,j)*geom::gs[j];
                tnsia[3*i+j]=geom::lu(i,j);
            }
        }
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            nsia[i]=spins::csx.getarrayelement(geom::lu(i,0),geom::lu(i,1),geom::lu(i,2));
        }
        CUDA_CALL(cudaMemcpy(Ccoord,tnsfa,3*geom::ss*sizeof(float),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Clu,tnsia,3*geom::ss*sizeof(int),cudaMemcpyHostToDevice));
        CUDA_CALL(cudaMemcpy(Czpsn,nsia,geom::ss*sizeof(int),cudaMemcpyHostToDevice));
        //zero the demag field array
        for(unsigned int i = 0 ; i < 3*geom::ss ; i++){tnsfa[i]=0.0;}CUDA_CALL(cudaMemcpy(CHDemag,tnsfa,3*geom::ss*sizeof(float),cudaMemcpyHostToDevice));
        //--------------------------------------------------------------------------------


        //make sure we clean up when the program exits
        atexit(deallocate_cuda_memory);
        cullbmi=true;
    }

    void setup_fourier_transform()
    {
        assert(intmat::imi);
        assert(intmat::imft);
        assert(cullbmi);

        //Even though we have 9 interaction matrices, 3 field arrays and
        //3 spin arrays we only need one transform in cufft. This is because
        //we can reuse the plan and alternate the sign depending on whether
        //we have a forward or a back transform
        /*Create a 3D FFT plan. */
        if(cufftPlan3d(&C3DP,geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],CUFFT_C2C)!=CUFFT_SUCCESS)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CUFFT 3D plan creation failed");
        }

        //At this point we can copy the interaction matrix from the CPU
        //as there is no need to do the determination of the interaction
        //matrix on the card. By asserting imft we are assuring that the
        //fourier transform of the interaction matrix has already been
        //completed.
        if(!intmat::imft)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You must perform the FT of the interaction matrix.");
        }
        else
        {

            //declare a holder on the heap
            Array3D<fftwf_complex> tempxx,tempxy,tempxz,tempyx,tempyy,tempyz,tempzx,tempzy,tempzz;
            tempxx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempxy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempxz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempyx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempyy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempyz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempzx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempzy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            tempzz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
            for(unsigned int i = 0 ; i < geom::zpdim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::zpdim[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::zpdim[2] ; k++)
                    {
                        for(unsigned int l = 0 ; l < 2 ; l++)
                        {
                            tempxx(i,j,k)[l]=float(intmat::Nxx(i,j,k)[l]);
                            tempxy(i,j,k)[l]=float(intmat::Nxy(i,j,k)[l]);
                            tempxz(i,j,k)[l]=float(intmat::Nxz(i,j,k)[l]);
                            tempyx(i,j,k)[l]=float(intmat::Nyx(i,j,k)[l]);
                            tempyy(i,j,k)[l]=float(intmat::Nyy(i,j,k)[l]);
                            tempyz(i,j,k)[l]=float(intmat::Nyz(i,j,k)[l]);
                            tempzx(i,j,k)[l]=float(intmat::Nzx(i,j,k)[l]);
                            tempzy(i,j,k)[l]=float(intmat::Nzy(i,j,k)[l]);
                            tempzz(i,j,k)[l]=float(intmat::Nzz(i,j,k)[l]);
                        }
                    }
                }
            }

            CUDA_CALL(cudaMemcpy(CCNxx,tempxx.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNxy,tempxy.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNxz,tempxz.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNyx,tempyx.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNyy,tempyy.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNyz,tempyz.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNzx,tempzx.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNzy,tempzy.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(CCNzz,tempzz.ptr(),geom::zps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
        }
        check_cuda_errors(__FILE__,__LINE__);
        cullbffti=true;
    }

    void deallocate_cuda_memory()
    {
        config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**EXIT information***" << std::endl;
        FIXOUT(config::Info,"Freeing space on GPU device" << std::flush);
        CUDA_CALL(cudaFree(CCNxx));
        CUDA_CALL(cudaFree(CCNxy));
        CUDA_CALL(cudaFree(CCNxz));
        CUDA_CALL(cudaFree(CCNyx));
        CUDA_CALL(cudaFree(CCNyy));
        CUDA_CALL(cudaFree(CCNyz));
        CUDA_CALL(cudaFree(CCNzx));
        CUDA_CALL(cudaFree(CCNzy));
        CUDA_CALL(cudaFree(CCNzz));
        CUDA_CALL(cudaFree(CCSx));
        CUDA_CALL(cudaFree(CCSy));
        CUDA_CALL(cudaFree(CCSz));
        CUDA_CALL(cudaFree(CCHx));
        CUDA_CALL(cudaFree(CCHy));
        CUDA_CALL(cudaFree(CCHz));
        CUDA_CALL(cudaFree(Cspin));
        CUDA_CALL(cudaFree(Cespin));
        CUDA_CALL(cudaFree(Cfspin));
        CUDA_CALL(cudaFree(CTemp));
        CUDA_CALL(cudaFree(Crand));
        CUDA_CALL(cudaFree(Cxadj));
        CUDA_CALL(cudaFree(Cadjncy));
        CUDA_CALL(cudaFree(Ccoord));
        CUDA_CALL(cudaFree(CHDemag));
        CUDA_CALL(cudaFree(Clu));
        CUDA_CALL(cudaFree(Czpsn));
		CUDA_CALL(cudaFree(CsurfArea));
		CUDA_CALL(cudaFree(CW1pf));
		CUDA_CALL(cudaFree(CW2pf));
		CUDA_CALL(cudaFree(Csigma));
        CUDA_CALL(cudaFree(Cfn));
        config::Info << "Done" << std::endl;
    }
}
