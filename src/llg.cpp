// File: llg.cpp
// Author:Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 21 Jan 2013 16:24:57
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
#include "../inc/mat.h"
namespace llg
{
    Array<double> fnx;
    Array<double> fny;
    Array<double> fnz;
    double applied[3]={0,0,0},T,dt,rdt,llgpf;
    void initLLG(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**LLG details***" << std::endl;
        FIXOUT(config::Info,"Initializing arrays:" << std::flush);
        fnx.resize(geom::nspins);
        fny.resize(geom::nspins);
        fnz.resize(geom::nspins);
        fnx.IFill(0);
        fny.IFill(0);
        fnz.IFill(0);

        SUCCESS(config::Info);
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

        libconfig::Setting &setting = config::cfg.lookup("llg");
        setting.lookupValue("dt",dt);
        for(unsigned int i = 0 ; i < 3 ;i++)
        {
            applied[i]=setting["applied"][i];
        }
        FIXOUTVEC(config::Info,"Applied field:",applied[0],applied[1],applied[2]);
        FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
        rdt=dt*mat::gamma;

        mat::sigma = sqrt(2.0*1.38e-23*mat::lambda/(mat::mu*mat::muB*dt*mat::gamma));
        FIXOUT(config::Info,"Sigma prefactor:" << mat::sigma << std::endl);
        FIXOUT(config::Info,"Resizing llg work arrays:" << std::flush);
        spins::eSx.resize(geom::nspins);
        spins::eSy.resize(geom::nspins);
        spins::eSz.resize(geom::nspins);
        fields::Hthx.resize(geom::nspins);
        fields::Hthy.resize(geom::nspins);
        fields::Hthz.resize(geom::nspins);
        SUCCESS(config::Info);
        llgpf = -1./(1.0+mat::lambda*mat::lambda);
    }

    void llgCPU(unsigned int t)
    {
        //calculate the 2 spin fields (dipolar, exchange, anisotropy)
        fields::ftdip();
        const double sqrtT=sqrt(T);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            fields::Hthx[i]=sqrtT*mat::sigma*Random::normal();
            fields::Hthy[i]=sqrtT*mat::sigma*Random::normal();
            fields::Hthz[i]=sqrtT*mat::sigma*Random::normal();

        }

        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            double h[3]={applied[0]+fields::Hthx[i]+fields::Hx[i],applied[1]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthx[i]+applied[2]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            fnx[i]=llgpf*(sxh[0]+mat::lambda*sxsxh[0]);
            fny[i]=llgpf*(sxh[1]+mat::lambda*sxsxh[1]);
            fnz[i]=llgpf*(sxh[2]+mat::lambda*sxsxh[2]);
            spins::eSx[i] = s[0] + fnx[i]*rdt;
            spins::eSy[i] = s[1] + fny[i]*rdt;
            spins::eSz[i] = s[2] + fnz[i]*rdt;
        }
        //perform the calculation of the 2 spin fields using the euler spin arrays
        fields::eftdip();
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={applied[0]+fields::Hthx[i]+fields::Hx[i],applied[1]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthx[i]+applied[2]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            const double fnp1[3]={llgpf*(sxh[0]+mat::lambda*sxsxh[0]),llgpf*(sxh[1]+mat::lambda*sxsxh[1]),llgpf*(sxh[2]+mat::lambda*sxsxh[2])};

            spins::Sx[i] +=(0.5*(fnx[i]+fnp1[0])*rdt);
            spins::Sy[i] +=(0.5*(fny[i]+fnp1[1])*rdt);
            spins::Sz[i] +=(0.5*(fnz[i]+fnp1[2])*rdt);

        }
    }
}
