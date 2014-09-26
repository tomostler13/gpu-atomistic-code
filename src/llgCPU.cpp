// File: llg.cpp
// Author:Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 26 Sep 2014 15:52:59
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include "../inc/arrays.h"
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
namespace llgCPU
{
    Array<double> fnx;
    Array<double> fny;
    Array<double> fnz;
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
        //calculate the 2 spin fields (dipolar, exchange)
        fields::ftdip();
        //FOR DEBUGGING THE FIELD
        /*        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                  {
                  std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
                  }
                  exit(0);*/
        const double sqrtT=sqrt(llg::T);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            fields::Hthx[i]=sqrtT*geom::sigma[i]*Random::normal();
            fields::Hthy[i]=sqrtT*geom::sigma[i]*Random::normal();
            fields::Hthz[i]=sqrtT*geom::sigma[i]*Random::normal();

        }
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            double h[3]={llg::applied[0]+fields::Hthx[i]+fields::Hx[i],llg::applied[1]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthx[i]+llg::applied[2]};
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
        //perform the calculation of the 2 spin fields using the euler spin arrays
        fields::eftdip();
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={llg::applied[0]+fields::Hthx[i]+fields::Hx[i],llg::applied[1]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthx[i]+llg::applied[2]};
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
