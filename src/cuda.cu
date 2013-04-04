// File: cuda.cu
// Author:Tom Ostler
// Last-modified: 04 Apr 2013 15:26:09
// Formally cuLLB.cu
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
//floating point version
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
	static  cufftComplex *CCNxx=NULL;
	static  cufftComplex *CCNxy=NULL;
	static  cufftComplex *CCNxz=NULL;
	static  cufftComplex *CCNyx=NULL;
	static  cufftComplex *CCNyy=NULL;
	static  cufftComplex *CCNyz=NULL;
	static  cufftComplex *CCNzx=NULL;
	static  cufftComplex *CCNzy=NULL;
	static  cufftComplex *CCNzz=NULL;
	static  cufftComplex *CCSkx=NULL;
	static  cufftComplex *CCSky=NULL;
	static  cufftComplex *CCSkz=NULL;
	static  cufftReal *CCSrx=NULL;
	static  cufftReal *CCSry=NULL;
	static  cufftReal *CCSrz=NULL;
	static  cufftComplex *CCHkx=NULL;
	static  cufftComplex *CCHky=NULL;
	static  cufftComplex *CCHkz=NULL;
	static  cufftReal *CCHrx=NULL;
	static  cufftReal *CCHry=NULL;
	static  cufftReal *CCHrz=NULL;

	//device pointers
	static  double *Cspin=NULL;
	static  double *Cespin=NULL;
	static  float *Crand=NULL;
	static  float *CH=NULL;
	static  int *Czpsn=NULL;//The is the zero pad spin number
	static  int *Clu=NULL;
	static  double *Cfn=NULL;
	//cufft plans
	cufftHandle C3DPr2c,C3DPc2r;

	void initGPU()
	{
//		CUDA_CALL(cudaDeviceReset());
	}
	void llgGPU(unsigned int& t)
	{

		//copy the spin data to the zero padded arrays
		cufields::CCopySpin<<<zpblockspergrid,threadsperblock>>>(geom::zps,geom::nspins,Cspin,Clu,CCSrx,CCSry,CCSrz,CCHrx,CCHry,CCHrz);
		//forward transform
		spins_forward();
		//perform convolution
		cufields::CFConv<<<czpblockspergrid,threadsperblock>>>(geom::czps,CCNxx,CCNxy,CCNxz,CCNyx,CCNyy,CCNyz,CCNzx,CCNzy,CCNzz,CCHkx,CCHky,CCHkz,CCSkx,CCSky,CCSkz);
		//transform the fields back
		fields_back();
		//copy the fields from the zero padded array to the demag field array
		cufields::CCopyFields<<<blockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CH,Czpsn,CCHrx,CCHry,CCHrz);
/*		Array3D<float> temp2x(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
		Array3D<float> temp2y(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
		Array3D<float> temp2z(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);;
		cudaMemcpy(temp2x.ptr(),CCHrx,geom::zps*sizeof(float),cudaMemcpyDeviceToHost);
		cudaMemcpy(temp2y.ptr(),CCHry,geom::zps*sizeof(float),cudaMemcpyDeviceToHost);
		cudaMemcpy(temp2z.ptr(),CCHrz,geom::zps*sizeof(float),cudaMemcpyDeviceToHost);*/
/*		//FOR DEBUGGING THE DIPOLAR FIELD/
		float temp1[3*geom::nspins];
		CUDA_CALL(cudaMemcpy(temp1,CH,3*geom::nspins*sizeof(float),cudaMemcpyDeviceToHost));
		for(unsigned int i = 0 ; i < geom::nspins ; i++)
		{
			int ijk[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
			std::cout << i << "\t" << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << "\t" << temp1[3*i] << "\t" << temp1[3*i+1] << "\t" << temp1[3*i+2] << std::endl;
			//std::cerr << temp2x.getarrayelement(ijk[0],ijk[1],ijk[2]) << std::endl;
			//std::cout << i << "\t" << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << "\t" << temp2x(ijk[0],ijk[1],ijk[2])/double(geom::zps) << "\t" << temp2y(ijk[0],ijk[1],ijk[2])/double(geom::zps) << "\t" << temp2z(ijk[0],ijk[1],ijk[2])/double(geom::zps) << std::endl;

		}
		exit(0);*/

		//generate the random numbers
		CURAND_CALL(curandGenerateNormal(gen,Crand,3*geom::nspins,0.0,1.0));
		cuint::CHeun1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,mat::sigma,llg::llgpf,mat::lambda,llg::rdt,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,Crand,Cfn);
		cufields::CCopySpin<<<zpblockspergrid,threadsperblock>>>(geom::zps,geom::nspins,Cespin,Clu,CCSrx,CCSry,CCSrz,CCHrx,CCHry,CCHrz);
		//forward transform
		spins_forward();
		//perform convolution
		cufields::CFConv<<<czpblockspergrid,threadsperblock>>>(geom::czps,CCNxx,CCNxy,CCNxz,CCNyx,CCNyy,CCNyz,CCNzx,CCNzy,CCNzz,CCHkx,CCHky,CCHkz,CCSkx,CCSky,CCSkz);
		//transform the fields back
		fields_back();
		//copy the fields from the zero padded array to the demag field array
		cufields::CCopyFields<<<blockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CH,Czpsn,CCHrx,CCHry,CCHrz);

		cuint::CHeun2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,mat::sigma,llg::llgpf,mat::lambda,llg::rdt,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,Crand,Cfn);
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
		/*        if(deviceProp.major < 2)
				  {
				  error::errPreamble(__FILE__,__LINE__);
				  error::errMessage("Cuda compute capability of 2.0 or greater is required for the use of function pointers.");
				  }
		 */
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
/*		if((cudaSetDevice(device))!=cudaSuccess)
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
		}*/

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
		Array3D<fftwf_complex> tempxx,tempxy,tempxz,tempyx,tempyy,tempyz,tempzx,tempzy,tempzz;
		tempxx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempxy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempxz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempyx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempyy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempyz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempzx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempzy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		tempzz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
		for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
		{
			for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
			{
				for(unsigned int k = 0 ; k < geom::cplxdim ; k++)
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

		CUDA_CALL(cudaMemcpy(CCNxx,tempxx.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNxy,tempxy.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNxz,tempxz.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNyx,tempyx.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNyy,tempyy.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNyz,tempyz.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNzx,tempzx.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNzy,tempzy.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(CCNzz,tempzz.ptr(),geom::czps*sizeof(cufftComplex),cudaMemcpyHostToDevice));
		intmat::Nxx.clear();
		intmat::Nxy.clear();
		intmat::Nxz.clear();
		intmat::Nyx.clear();
		intmat::Nyy.clear();
		intmat::Nyz.clear();
		intmat::Nzx.clear();
		intmat::Nzy.clear();
		intmat::Nzz.clear();
		tempxx.clear();
		tempxy.clear();
		tempxz.clear();
		tempyx.clear();
		tempyy.clear();
		tempyz.clear();
		tempzx.clear();
		tempzy.clear();
		tempzz.clear();
		check_cuda_errors(__FILE__,__LINE__);
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
}
