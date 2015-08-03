// File: llg.cpp
// Author:Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 03 Aug 2015 17:00:18
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include "../inc/arrays.h"
#include "../inc/exch.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/intmat.h"
#include "../inc/defines.h"
#include "../inc/llg.h"
#include "../inc/random.h"
#include "../inc/anis.h"
#include "../inc/matrix_mul_cpu.h"
namespace llgCPU
{
    Array<double> fnx;
    Array<double> fny;
    Array<double> fnz;
    Array<double> Ts,cps,dps;
    void initLLG(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**LLG CPU details***" << std::endl;
        FIXOUT(config::Info,"Initializing arrays:" << std::flush);
        fnx.resize(geom::nspins);
        fny.resize(geom::nspins);
        fnz.resize(geom::nspins);
        fnx.IFill(0);
        fny.IFill(0);
        fnz.IFill(0);
        std::cout << "Resizing to " << geom::ucm.GetNMS() << std::endl;

        SUCCESS(config::Info);
        FIXOUT(config::Info,"Resizing llg (CPU) work arrays:" << std::flush);
        spins::eSx.resize(geom::nspins);
        spins::eSy.resize(geom::nspins);
        spins::eSz.resize(geom::nspins);
        fields::Hthx.resize(geom::nspins);
        fields::Hthy.resize(geom::nspins);
        fields::Hthz.resize(geom::nspins);
        SUCCESS(config::Info);
    }

    void llgCPU(unsigned int t)
    {
        fields::Hx.IFill(0);
        fields::Hy.IFill(0);
        fields::Hz.IFill(0);
        //calculate the 2 spin fields (dipolar, exchange)

        if(config::exchm==0 || config::exchm>98)
        {
            fields::ftdip();
        }
        else if(config::dipm==0 && config::inc_dip==true)
        {
            if(t%spins::update==0)
            {
                fields::dipftdip();
            }
        }
//        fields::bfdip();
//        else if(config::dipm==
        //FOR DEBUGGING THE FIELD
/*        if(t==0)
        {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::HDemagx(i) << "\t" << fields::HDemagy(i) << "\t" << fields::HDemagz(i) << std::endl;
        }
        exit(0);
        }*/
        //exit(0);*/
        //then we call the DIA SpMV multiplication
        if(config::exchm==1)
        {
            matmul::spmv_dia_diag(geom::nspins,exch::diagnumdiag,exch::diagoffset,exch::dataxx,exch::datayy,exch::datazz,spins::Sx,spins::Sy,spins::Sz,fields::Hx,fields::Hy,fields::Hz);
            if(config::offdiag==true)
            {
                matmul::spmv_dia_offdiag(geom::nspins,exch::offdiagnumdiag,exch::offdiagoffset,exch::dataxy,exch::dataxz,exch::datayx,exch::datayz,exch::datazx,exch::datazy,spins::Sx,spins::Sy,spins::Sz,fields::Hx,fields::Hy,fields::Hz);
            }
        }
        else if(config::exchm==2)
        {
            matmul::spmv_csr_diag(geom::nspins,exch::xadj,exch::adjncy,exch::dataxx,exch::datayy,exch::datazz,spins::Sx,spins::Sy,spins::Sz,fields::Hx,fields::Hy,fields::Hz);
            if(config::offdiag==true)
            {
                matmul::spmv_csr_offdiag();
            }
        }
        else if(config::exchm>98)
        {
            matmul::spmv_csr_diag(geom::nspins,exch::xadj,exch::adjncy,exch::dataxx,exch::datayy,exch::datazz,spins::Sx,spins::Sy,spins::Sz,fields::Hx,fields::Hy,fields::Hz);
            if(config::offdiag==true)
            {
                matmul::spmv_csr_offdiag();
            }
        }
        //FOR DEBUGGING THE FIELD
        /*if(t==0)
        {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
        }
        exit(0);
        }
        //exit(0);*/
        const double sqrtT=sqrt(llg::T);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            fields::Hthx[i]=sqrtT*geom::sigma[i]*Random::normal();
            fields::Hthy[i]=sqrtT*geom::sigma[i]*Random::normal();
            fields::Hthz[i]=sqrtT*geom::sigma[i]*Random::normal();
            //std::cout << fields::Hthx[i] << "\t" << fields::Hthy[i] << "\t" << fields::Hthz[i] << std::endl;
        }


        //calculate the four-spin exchange
        if(exch::inc4spin)
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                fields::H4sx[i]=0.0;
                fields::H4sy[i]=0.0;
                fields::H4sz[i]=0.0;
                double h[3]={0,0,0};
                for(unsigned int q = exch::xadj_j[i] ; q < exch::xadj_j[i+1] ; q++)
                {
                    unsigned int j=exch::adjncy_j[q];
                    unsigned int k=exch::adjncy_k[q];
                    unsigned int l=exch::adjncy_l[q];

                    //std::cout << "spin i=" << i << "\thas quartet\t" << j << "\t" << k << "\t" << l << std::endl;
                    const double sj[3]={spins::Sx[j],spins::Sy[j],spins::Sz[j]};
                    const double sk[3]={spins::Sx[k],spins::Sy[k],spins::Sz[k]};
                    const double sl[3]={spins::Sx[l],spins::Sy[l],spins::Sz[l]};
                    const double skdotsl=sk[0]*sl[0]+sk[1]*sl[1]+sk[2]*sl[2];
                    const double sjdotsl=sj[0]*sl[0]+sj[1]*sl[1]+sj[2]*sl[2];
                    const double skdotsj=sk[0]*sj[0]+sk[1]*sj[1]+sk[2]*sj[2];
                    fields::H4sx[i]-=(exch::JQ(0)*(sj[0]*skdotsl+sk[0]*sjdotsl+sl[0]*skdotsj));
                    fields::H4sy[i]-=(exch::JQ(0)*(sj[1]*skdotsl+sk[1]*sjdotsl+sl[1]*skdotsj));
                    fields::H4sz[i]-=(exch::JQ(0)*(sj[2]*skdotsl+sk[2]*sjdotsl+sl[2]*skdotsj));
//h[2]-=(exch::JQ(0,0)*(sj[2]*skdotsl+sk[2]*sjdotsl+sl[2]*skdotsj));
                }
//std::cout << h[0] << "\t" << h[1] << "\t" << h[2] << std::endl;
//std::cin.get();

                //std::cout << fields::H4sx[i] << "\t" << fields::H4sy[i] << "\t" << fields::H4sz[i] << std::endl;
            }
        }
//     exit(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            //Adding the demag here should be zero (see fields.cpp) if config::inc_dip==false
            double h[3]={llg::applied[0]+fields::Hthx[i]+fields::Hx[i]+fields::HDemagx[i]+fields::H4sx[i],llg::applied[1]+fields::Hthy[i]+fields::Hy[i]+fields::HDemagy[i]+fields::H4sy[i],fields::Hz[i]+fields::Hthx[i]+llg::applied[2]+fields::HDemagz[i]+fields::H4sz[i]};
/*        if(t==0)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << h[0] << "\t" << h[1] << "\t" << h[2] << std::endl;
        }*/
            //std::cout << h[0] << "\t" << h[1] << "\t" << h[2] << std::endl;
            //--------------------------------------------------------
            //calculate the first order uniaxial anisotropy component
            //direction of the axis
            const double d[3]={anis::k1udir(i,0),anis::k1udir(i,1),anis::k1udir(i,2)};
            // S . n (n=direction of anisotropy axis)
            const double sdn=s[0]*d[0]+s[1]*d[1]+s[2]*d[2];
            // 2 * D (where D is anisotropy constant
            const double TwoD=anis::k1u[i];
            h[0]+=(TwoD*sdn*d[0]);
            h[1]+=(TwoD*sdn*d[1]);
            h[2]+=(TwoD*sdn*d[2]);
            //--------------------------------------------------------

            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            fnx[i]=geom::llgpf[i]*(sxh[0]+geom::lambda[i]*sxsxh[0]);
            fny[i]=geom::llgpf[i]*(sxh[1]+geom::lambda[i]*sxsxh[1]);
            fnz[i]=geom::llgpf[i]*(sxh[2]+geom::lambda[i]*sxsxh[2]);
            spins::eSx[i] = s[0] + fnx[i]*llg::rdt;
            spins::eSy[i] = s[1] + fny[i]*llg::rdt;
            spins::eSz[i] = s[2] + fnz[i]*llg::rdt;
            const double mods = sqrt(spins::eSx[i]*spins::eSx[i]+spins::eSy[i]*spins::eSy[i]+spins::eSz[i]*spins::eSz[i]);
            spins::eSx[i]/=mods;
            spins::eSy[i]/=mods;
            spins::eSz[i]/=mods;
        }
        //exit(0);

        //perform the calculation of the 2 spin fields using the euler spin arrays
        if(config::exchm==0)
        {
            fields::eftdip();
        }
/*        else if(config::dipm==0 && config::inc_dip==true)
        {
            fields::dipeftdip();
        }*/
        //then we call the DIA SpMV multiplication
        if(config::exchm==1)
        {
            matmul::spmv_dia_diag(geom::nspins,exch::diagnumdiag,exch::diagoffset,exch::dataxx,exch::datayy,exch::datazz,spins::eSx,spins::eSy,spins::eSz,fields::Hx,fields::Hy,fields::Hz);
            if(config::offdiag==true)
            {
                matmul::spmv_dia_offdiag(geom::nspins,exch::offdiagnumdiag,exch::offdiagoffset,exch::dataxy,exch::dataxz,exch::datayx,exch::datayz,exch::datazx,exch::datazy,spins::eSx,spins::eSy,spins::eSz,fields::Hx,fields::Hy,fields::Hz);
            }
        }
        else if(config::exchm==2)
        {
            matmul::spmv_csr_diag(geom::nspins,exch::xadj,exch::adjncy,exch::dataxx,exch::datayy,exch::datazz,spins::eSx,spins::eSy,spins::eSz,fields::Hx,fields::Hy,fields::Hz);
            if(config::offdiag==true)
            {
                matmul::spmv_csr_offdiag();
            }
        }
        if(exch::inc4spin)
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                fields::H4sx[i]=0.0;
                fields::H4sy[i]=0.0;
                fields::H4sz[i]=0.0;
                for(unsigned int q = exch::xadj_j[i] ; q < exch::xadj_j[i+1] ; q++)
                {
                    unsigned int j=exch::adjncy_j[q];
                    unsigned int k=exch::adjncy_k[q];
                    unsigned int l=exch::adjncy_l[q];
                    const double sj[3]={spins::eSx[j],spins::eSy[j],spins::eSz[j]};
                    const double sk[3]={spins::eSx[k],spins::eSy[k],spins::eSz[k]};
                    const double sl[3]={spins::eSx[l],spins::eSy[l],spins::eSz[l]};
                    const double skdotsl=sk[0]*sl[0]+sk[1]*sl[1]+sk[2]*sl[2];
                    const double sjdotsl=sj[0]*sl[0]+sj[1]*sl[1]+sj[2]*sl[2];
                    const double skdotsj=sk[0]*sj[0]+sk[1]*sj[1]+sk[2]*sj[2];
                    fields::H4sx[i]-=(exch::JQ(0)*(sj[0]*skdotsl+sk[0]*sjdotsl+sl[0]*skdotsj));
                    fields::H4sy[i]-=(exch::JQ(0)*(sj[1]*skdotsl+sk[1]*sjdotsl+sl[1]*skdotsj));
                    fields::H4sz[i]-=(exch::JQ(0)*(sj[2]*skdotsl+sk[2]*sjdotsl+sl[2]*skdotsj));
                }
            }
        }
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={llg::applied[0]+fields::Hthx[i]+fields::Hx[i]+fields::HDemagx[i]+fields::H4sx[i],llg::applied[1]+fields::Hthy[i]+fields::Hy[i]+fields::HDemagy[i]+fields::H4sy[i],fields::Hz[i]+fields::Hthx[i]+llg::applied[2]+fields::HDemagz[i]+fields::H4sz[i]};
            //--------------------------------------------------------
            //calculate the first order uniaxial anisotropy component
            //direction of the axis
            const double d[3]={anis::k1udir(i,0),anis::k1udir(i,1),anis::k1udir(i,2)};
            // S . n (n=direction of anisotropy axis)
            const double sdn=s[0]*d[0]+s[1]*d[1]+s[2]*d[2];
            // 2 * D (where D is anisotropy constant
            const double TwoD=anis::k1u[i];
            h[0]+=(TwoD*sdn*d[0]);
            h[1]+=(TwoD*sdn*d[1]);
            h[2]+=(TwoD*sdn*d[2]);
            //--------------------------------------------------------

            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            const double fnp1[3]={geom::llgpf[i]*(sxh[0]+geom::lambda[i]*sxsxh[0]),geom::llgpf[i]*(sxh[1]+geom::lambda[i]*sxsxh[1]),geom::llgpf[i]*(sxh[2]+geom::lambda[i]*sxsxh[2])};

            spins::Sx[i] +=(0.5*(fnx[i]+fnp1[0])*llg::rdt);
            spins::Sy[i] +=(0.5*(fny[i]+fnp1[1])*llg::rdt);
            spins::Sz[i] +=(0.5*(fnz[i]+fnp1[2])*llg::rdt);
            const double mods = sqrt(spins::Sx[i]*spins::Sx[i]+spins::Sy[i]*spins::Sy[i]+spins::Sz[i]*spins::Sz[i]);
            spins::Sx[i]/=mods;
            spins::Sy[i]/=mods;
            spins::Sz[i]/=mods;
        }

    }
}
