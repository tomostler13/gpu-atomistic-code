// File: cuintRK5.cu
// Author:Tom Ostler
// Last-modified: 19 May 2023 14:10:44
#include "../inc/cufields.h"
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/util.h"
#include "../inc/cuintRK5.h"
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
namespace cuintRK5
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

    // Constants from Romeo et al. Physica B Vol. 403 464-468 (2008)
    __constant__ double Ai[6]={0.0,0.2,0.3,0.6,1.0,0.875};
    __constant__ double Ci[6]={37.0/378.0,0.0,250.0/621.0,125.0/594.0,0.0,512.0/1771.0};
    __constant__ double Cpi[6]={2825.0/27648.0,0.0,18575.0/48384.0,13525.0/55296.0,277.0/14336.0,0.25};
    __constant__ double B21 = 0.2;
    __constant__ double B3j[2]={3.0/40.0,9.0/40.0};
    __constant__ double B4j[3]={0.3,-0.9,1.2};
    __constant__ double B5j[4]={-11.0/54,2.5,-70.0/27.0,35.0/27.0};
    __constant__ double B6j[5]={1631.0/55296.0,175.0/512.0,575.0/13824.0,44275.0/110592.0,253.0/4096.0};


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
    __global__ void CdetRK5k1(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRKk1,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
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

            double h[3]={CH[3*i]+CHstg[3*i]+lrn[0]+appliedx,CH[3*i+1]+CHstg[3*i+1]+lrn[1]+appliedy,CH[3*i+2]+CHstg[3*i+2]+lrn[2]+appliedz};

            //printf("%4.5f\t%4.5f\t%4.5f\n",CH[3*i],CH[3*i+1],CH[3*i+2]);

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
            double mods=0.0;
            double es[3]={0.0,0.0,0.0};
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lk1[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
                es[j]=s[j]+Crdt*(B21*lk1[j]);
                mods+=es[j]*es[j];
            }
            //calculate one over the square root of the spin modulus
            const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //set the euler spin value
                Cespin[3*i+j]=es[j]*nf;
                CRKk1[3*i+j]=lk1[j];
            }
        }
    }

    __global__ void CdetRK5k2(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRKk2,double *CRKk1,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
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

            double h[3]={CH[3*i]+CHstg[3*i]+lrn[0]+appliedx,CH[3*i+1]+CHstg[3*i+1]+lrn[1]+appliedy,CH[3*i+2]+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            double es[3]={Cespin[3*i],Cspin[3*i+1],Cespin[3*i+2]};
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
            double mods=0.0;
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lk2[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
                //Cfn[3*i+j]=lfn[j];
                es[j]=s[j]+Crdt*(B3j[0]*CRKk1[3*i+j]+B3j[1]*lk2[j]);
                mods+=es[j]*es[j];
            }
            //calculate one over the square root of the spin modulus
            const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //set the euler spin value
                Cespin[3*i+j]=es[j]*nf;
                CRKk2[3*i+j]=lk2[j];
            }
        }
    }

    __global__ void CdetRK5k3(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRKk3,double *CRKk2,double *CRKk1,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
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

            double h[3]={CH[3*i]+CHstg[3*i]+lrn[0]+appliedx,CH[3*i+1]+CHstg[3*i+1]+lrn[1]+appliedy,CH[3*i+2]+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            double es[3]={Cespin[3*i],Cspin[3*i+1],Cespin[3*i+2]};
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

            double lk3[3]={0,0,0};
            double mods=0.0;
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lk3[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
                //Cfn[3*i+j]=lfn[j];
                es[j]=s[j] + Crdt*(B4j[0]*CRKk1[3*i+j] + B4j[1]*CRKk2[3*i+j] + B4j[2]*lk3[j]);
                mods+=es[j]*es[j];
            }
            //calculate one over the square root of the spin modulus
            const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //set the euler spin value
                Cespin[3*i+j]=es[j]*nf;
                CRKk3[3*i+j]=lk3[j];
            }
        }
    }
    __global__ void CdetRK5k4(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRKk4,double *CRKk3,double *CRKk2,double *CRKk1,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
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

            double h[3]={CH[3*i]+CHstg[3*i]+lrn[0]+appliedx,CH[3*i+1]+CHstg[3*i+1]+lrn[1]+appliedy,CH[3*i+2]+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            double es[3]={Cespin[3*i],Cspin[3*i+1],Cespin[3*i+2]};
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

            double lk4[3]={0,0,0};
            double mods=0.0;
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lk4[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
                es[j]=s[j]+Crdt*(B5j[0]*CRKk1[3*i+j] + B5j[1]*CRKk2[3*i+j] + B5j[2]*CRKk3[3*i+j] + B5j[3]*lk4[j]);
                mods+=es[j]*es[j];
            }
            //calculate one over the square root of the spin modulus
            const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //set the euler spin value
                Cespin[3*i+j]=es[j]*nf;
                CRKk4[3*i+j]=lk4[j];
            }
        }
    }
    __global__ void CdetRK5k5(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRKk5,double *CRKk4,double *CRKk3,double *CRKk2,double *CRKk1,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
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

            double h[3]={CH[3*i]+CHstg[3*i]+lrn[0]+appliedx,CH[3*i+1]+CHstg[3*i+1]+lrn[1]+appliedy,CH[3*i+2]+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            double es[3]={Cespin[3*i],Cspin[3*i+1],Cespin[3*i+2]};
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

            double lk5[3]={0,0,0};
            double mods=0.0;
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lk5[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
                es[j]=s[j]+Crdt*(B6j[0]*CRKk1[3*i+j] + B6j[1]*CRKk2[3*i+j] + B6j[2]*CRKk3[3*i+j] + B6j[3]*CRKk4[3*i+j] + B6j[4]*lk5[j]);
                mods+=es[j]*es[j];
            }
            //calculate one over the square root of the spin modulus
            const double nf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                //set the euler spin value
                Cespin[3*i+j]=es[j]*nf;
                CRKk5[3*i+j]=lk5[j];
            }
        }
    }
    __global__ void CdetRK5k6(int N,double T,double appliedx,double appliedy,double appliedz,double *CH,double *Cspin,double *Cespin,double *CRKk6,double *CRKk5,double *CRKk4,double *CRKk3,double *CRKk2,double *CRKk1,float *Crand,double *Csigma,double *Cllgpf,double *Clambda,double *Ck1u,double *Ck1udir,double *CHstg)
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

            double h[3]={CH[3*i]+CHstg[3*i]+lrn[0]+appliedx,CH[3*i+1]+CHstg[3*i+1]+lrn[1]+appliedy,CH[3*i+2]+CHstg[3*i+2]+lrn[2]+appliedz};
            //printf("%4.5f\t%4.5f\t%4.5f\n",CHstg[3*i],CHstg[3*i+1],CHstg[3*i+2]);

            const double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            double es[3]={Cespin[3*i],Cspin[3*i+1],Cespin[3*i+2]};
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

            double lk6[3]={0,0,0};
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                lk6[j] = llgpf*(sxh[j]+lambda*sxsxh[j]);
            }
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                CRKk6[3*i+j]=lk6[j];
            }
        }
    }

    __global__ void CdetSnp1(int N,double *Cspin,double *CRKk1,double *CRKk2,double *CRKk3,double *CRKk4,double *CRKk5,double *CRKk6)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            double mods=0.0;
            double ls[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                ls[j]=ls[j]+ Crdt*(CRKk1[3*i+j]*Ci[0] + CRKk2[3*i+j]*Ci[1] + CRKk3[3*i+j]*Ci[2] + CRKk4[3*i+j]*Ci[3]+CRKk5[3*i+j]*Ci[4]+CRKk6[3*i+j]*Ci[5]);
                mods+=ls[j]*ls[j];
            }
            const double nrf=rsqrt(mods);
            for(unsigned int j = 0 ; j < 3 ; j++)
            {

//                printf("Spin %d Component %d Cspin %4.10f k1 %4.10f k2 %4.10f k3 %4.10f k4 %4.10f \n",i,j,Cspin[3*i+j],CRKk1[3*i+j],CRKk2[3*i+j],CRKk3[3*i+j],CRKk4[3*i+j]);
                Cspin[3*i+j] = ls[j]*nrf;
            }
        }
    }
}
