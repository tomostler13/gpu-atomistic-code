// File: cuda.cu
// Author:Tom Ostler
// Last-modified: 18 May 2023 05:09:16 PM
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
#include "../inc/cuintRK4.h"
#include "../inc/cuint.h"
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


    /*void rampStag(double rampval)
    {
        cufields::CincStagFields<<<blockspergrid,threadsperblock>>>(geom::nspins,rampval,CInitHstg,CHstg);
    }*/

    void llgGPURK4(unsigned int& t)
    {

        CURAND_CALL(curandGenerateNormal(gen,Crand,curandN,0.0,1.0));
        //in this case we are using the interaction matrix for both the exchange and the
        //dipole-dipole field so we might aswell update both at once
        if(config::exchm==0)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("GPU code for RK4 integration with FFT calculated exchange fields not yet implemented");
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
                //cufields::CdipCopyFields<<<zpblockspergrid,threadsperblock>>>(geom::nspins,geom::zps,CHDemag,CHr,Ckx,Cky,Ckz);
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
            //Calculate the fields based on the current spin value
            if(config::exchm==1)//DIA
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cspin,CHDemag,CH);
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);

                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                }
            }
            else if(config::exchm==2)//CSR
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cspin,CHDemag,CH);
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
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
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_DIA_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cdiagoffset,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                }
            }
            else if(config::exchm==2)//CSR
            {
                if(!config::offdiag)
                {
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cspin,CHDemag,CH);
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdyy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                }
                else if(config::offdiag)
                {
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cspin,CHDemag,CH);
                    cuintRK4::CdetRK4k1<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k1,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    //Here, we overwrite CH with a new field updated based on the S_n+dt*k1/2 field
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k2<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k2,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k3<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k3,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                    cufields::CSpMV_CSR_OffDiag<<<blockspergrid,threadsperblock>>>(geom::nspins,Cxadj,Cadjncy,Cdxx,Cdxy,Cdxz,Cdyx,Cdyy,Cdyz,Cdzx,Cdzy,Cdzz,Cespin,CHDemag,CH);
                    cuintRK4::CdetRK4k4<<<blockspergrid,threadsperblock>>>(geom::nspins,llg::T,llg::applied[0],llg::applied[1],llg::applied[2],CH,Cspin,Cespin,CRK4k4,Crand,Csigma,Cllgpf,Clambda,Ck1u,Ck1udir,CHstg);
                }
            }
        }
        cuintRK4::CdetSnp1<<<blockspergrid,threadsperblock>>>(geom::nspins,Cspin,CRK4k1,CRK4k2,CRK4k3,CRK4k4);
        //generate the random numbers
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
        //in this case we are using the interaction matrix for both the exchange and the
        //dipole-dipole field so we might aswell update both at once

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

}
