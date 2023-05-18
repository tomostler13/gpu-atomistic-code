// File: cuintRK4.cu
// Author:Tom Ostler
// Last-modified: 18 May 2023 03:27:02 PM
#include "../inc/cufields.h"
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/util.h"
#include "../inc/cuintRK4.h"
#include "../inc/fields.h"
#include "../inc/defines.h"
#include "../inc/cudadefs.h"
#include "../inc/anis.h"
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
namespace cuintRK4
{
    //Reduced timestep
    __constant__ double Crdt;
    __constant__ double CK2perp;
    __constant__ double CK2perpdir[3];
    __constant__ double CK2par;
    __constant__ double CK2pardir[3];
    __constant__ double CK4perp;
    __constant__ double CK4perpdirs[3][3]; 
    __constant__ double CK4par;
    __constant__ double CK4pardirs[3][3];
    __constant__ double Cbasemm;
    __constant__ double onesixth=1.0/6.0;
    __constant__ double onethird=1.0/3.0;


    void copyConstData()
    {
        FIXOUT(config::Info,"Copying const data with cuint scope to card:" << std::flush);
        cudaMemcpyToSymbol(*(&Crdt),&llg::rdt,sizeof(double));
        cudaMemcpyToSymbol(*(&CK2par),&anis::k2par,sizeof(double));
        cudaMemcpyToSymbol(*(&CK2perp),&anis::k2perp,sizeof(double));
        cudaMemcpyToSymbol(*(&CK4par),&anis::k4par,sizeof(double));
        cudaMemcpyToSymbol(*(&CK4perp),&anis::k4perp,sizeof(double));
        cudaMemcpyToSymbol(CK2pardir,anis::k2pardir.ptr(),3*sizeof(double));
        cudaMemcpyToSymbol(CK2perpdir,anis::k2perpdir.ptr(),3*sizeof(double));
        cudaMemcpyToSymbol(CK4pardirs,anis::k4pardirs.ptr(),9*sizeof(double));
        cudaMemcpyToSymbol(CK4perpdirs,anis::k4perpdirs.ptr(),9*sizeof(double));
        double bmm=geom::ucm.GetMu(0)*llg::muB;
        cudaMemcpyToSymbol(*(&Cbasemm),&bmm,sizeof(double));
        config::Info << "Done" << std::endl;
    }
    __global__ void CdetRK4k1(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRK4k1,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //The prefactor for the thermal term
            const double TP=sqrt(T)*Csigma[i];
            const double llgpf=Cllgpf[i];
            const double lambda=Clambda[i];
            const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
            //printf("%4.10f\t%4.10f\t%4.10f\n",lrn[0],lrn[1],lrn[2]);

            double h[3]={double(CH[3*i])+CHstg[3*i]+lrn[0]+appliedx,double(CH[3*i+1])+CHstg[3*i+1]+lrn[1]+appliedy,double(CH[3*i+2])+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            //calculate the field arising from the first order uniaxial anisotropy
            const double k1udir[3]={Ck1udir[3*i],Ck1udir[3*i+1],Ck1udir[3*i+2]};
            const double k1u=Ck1u[i];
            const double sdn = s[0]*k1udir[0] + s[1]*k1udir[1] + s[2]*k1udir[2];

            h[0]+=(k1u*sdn*k1udir[0]);
            h[1]+=(k1u*sdn*k1udir[1]);
            h[2]+=(k1u*sdn*k1udir[2]);

            //field from k2_perp
            const double sdk2perpdir=s[0]*CK2perpdir[0]+s[1]*CK2perpdir[1]+s[2]*CK2perpdir[2];
            h[0]+=(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[0];
            h[1]+=(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[1];
            h[2]+=(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[2];
            /*if(i==0)
            {
                printf("CK2perp %4.5f\t%4.5f\t%4.5f\n",(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[0],(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[1],(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[2]);
            }*/
            //field from k2_parallel
            const double sdk2pardir=s[0]*CK2pardir[0]+s[1]*CK2pardir[1]+s[2]*CK2pardir[2];
            h[0]+=(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[0];
            h[1]+=(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[1];
            h[2]+=(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[2];
            /*if(i==0)
            {
                printf("CK2par %4.5f\t%4.5f\t%4.5f\n",(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[0],(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[1],(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[2]);
            }*/

            //field from k4 perp and par
            for(unsigned int dir = 0 ; dir < 3; dir++)
            {
                const double sdk4perpdir=s[0]*CK4perpdirs[dir][0]+s[1]*CK4perpdirs[dir][1]+s[2]*CK4perpdirs[dir][2];
                const double sdk4pardir=s[0]*CK4pardirs[dir][0]+s[1]*CK4pardirs[dir][1]+s[2]*CK4pardirs[dir][2];
                h[0]+=(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][0];
                h[1]+=(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][1];
                h[2]+=(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][2];
                /*if(i==0)
                {
                    printf("CK4perp %4.5f\t%4.5f\t%4.5f\n",(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][0],(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][1],(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][2]);
                }*/


                h[0]+=(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][0];
                h[1]+=(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][1];
                h[2]+=(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][2];
                /*if(i==0)
                {
                    printf("CK4par %4.5f\t%4.5f\t%4.5f\n",(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][0],(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][1],(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][2]);
                }*/
            }


            //printf("%4.10f\t%4.10f\t%4.10f\n",h[0],h[1],h[2]);

            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            double lk1[3]={0,0,0};
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lk1[j] = Crdt*llgpf*(sxh[j]+lambda*sxsxh[j]);
                //Cfn[3*i+j]=lfn[j];
                //es[j]=s[j]+lfn[j]*Crdt;
                //mods+=es[j]*es[j];
            }
            //calculate one over the square root of the spin modulus
            //const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //set the euler spin value
                Cespin[3*i+j]=s[j]+Crdt*lk1[j];
                CRK4k1[3*i+j]=lk1[j];
            }
        }
    }
    __global__ void CdetRK4kn(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRK4k2,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            //The prefactor for the thermal term
            const double TP=sqrt(T)*Csigma[i];
            const double llgpf=Cllgpf[i];
            const double lambda=Clambda[i];
            const double lrn[3]={double(Crand[3*i])*TP,double(Crand[3*i+1])*TP,double(Crand[3*i+2])*TP};
            //printf("%4.10f\t%4.10f\t%4.10f\n",lrn[0],lrn[1],lrn[2]);

            double h[3]={double(CH[3*i])+CHstg[3*i]+lrn[0]+appliedx,double(CH[3*i+1])+CHstg[3*i+1]+lrn[1]+appliedy,double(CH[3*i+2])+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            const double es[3]={Cespin[3*i],Cspin[3*i+1],Cespin[3*i+2]};
            //calculate the field arising from the first order uniaxial anisotropy
            const double k1udir[3]={Ck1udir[3*i],Ck1udir[3*i+1],Ck1udir[3*i+2]};
            const double k1u=Ck1u[i];
            const double sdn = es[0]*k1udir[0] + es[1]*k1udir[1] + es[2]*k1udir[2];

            h[0]+=(k1u*sdn*k1udir[0]);
            h[1]+=(k1u*sdn*k1udir[1]);
            h[2]+=(k1u*sdn*k1udir[2]);

            //field from k2_perp
            const double sdk2perpdir=es[0]*CK2perpdir[0]+es[1]*CK2perpdir[1]+es[2]*CK2perpdir[2];
            h[0]+=(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[0];
            h[1]+=(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[1];
            h[2]+=(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[2];
            /*if(i==0)
            {
                printf("CK2perp %4.5f\t%4.5f\t%4.5f\n",(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[0],(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[1],(2.0*CK2perp/Cbasemm)*sdk2perpdir*CK2perpdir[2]);
            }*/
            //field from k2_parallel
            const double sdk2pardir=es[0]*CK2pardir[0]+es[1]*CK2pardir[1]+es[2]*CK2pardir[2];
            h[0]+=(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[0];
            h[1]+=(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[1];
            h[2]+=(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[2];
            /*if(i==0)
            {
                printf("CK2par %4.5f\t%4.5f\t%4.5f\n",(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[0],(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[1],(2.0*CK2par/Cbasemm)*sdk2pardir*CK2pardir[2]);
            }*/

            //field from k4 perp and par
            for(unsigned int dir = 0 ; dir < 3; dir++)
            {
                const double sdk4perpdir=es[0]*CK4perpdirs[dir][0]+es[1]*CK4perpdirs[dir][1]+es[2]*CK4perpdirs[dir][2];
                const double sdk4pardir=es[0]*CK4pardirs[dir][0]+es[1]*CK4pardirs[dir][1]+es[2]*CK4pardirs[dir][2];
                h[0]+=(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][0];
                h[1]+=(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][1];
                h[2]+=(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][2];
                /*if(i==0)
                {
                    printf("CK4perp %4.5f\t%4.5f\t%4.5f\n",(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][0],(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][1],(2.0*CK4perp/Cbasemm)*sdk4perpdir*CK4perpdirs[dir][2]);
                }*/


                h[0]+=(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][0];
                h[1]+=(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][1];
                h[2]+=(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][2];
                /*if(i==0)
                {
                    printf("CK4par %4.5f\t%4.5f\t%4.5f\n",(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][0],(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][1],(2.0*CK4par/Cbasemm)*sdk4pardir*CK4pardirs[dir][2]);
                }*/
            }


            //printf("%4.10f\t%4.10f\t%4.10f\n",h[0],h[1],h[2]);

            const double sxh[3]={es[1]*h[2] - es[2]*h[1],es[2]*h[0]-es[0]*h[2],es[0]*h[1]-es[1]*h[0]};
            const double sxsxh[3]={es[1]*sxh[2]-es[2]*sxh[1],es[2]*sxh[0]-es[0]*sxh[2],es[0]*sxh[1]-es[1]*sxh[0]};

            double lk2[3]={0,0,0};
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //This should be k2 from the RK4 integration...
                lk2[j] = Crdt*llgpf*(sxh[j]+lambda*sxsxh[j]);
                //Cfn[3*i+j]=lfn[j];
                //es[j]=s[j]+lfn[j]*Crdt;
                //mods+=es[j]*es[j];
            }
            //calculate one over the square root of the spin modulus
            //const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //set the euler spin value
                Cespin[3*i+j]=s[j]+Crdt*lk2[j];
                CRK4k2[3*i+j]=lk2[j];
            }
        }
    }

    __global__ void CdetSnp1(int N,double *Cspin,double *CRK4k1,double *CRK4k2,double *CRK4k3,double *CRK4k4)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                Cspin[3*i+j] = Cspin[3*i+j] + CRK4k1[3*i+j]*onesixth + CRK4k2[3*i+j]*onethird + CRK4k3[3*i+j]*onethird + CRK4k4[3*i+j]*onesixth;
            }
        }
    }
}
