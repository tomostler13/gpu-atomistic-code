// File: stepchange.cpp
// Author: Tom Ostler // Created: 29 Mar 2013
// Last-modified: 18 Oct 2016 13:15:13
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "../inc/arrays.h"
#include "../inc/defines.h"
#include "../inc/util.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/intmat.h"
#include "../inc/fields.h"
#include "../inc/llg.h"
#include "../inc/sim.h"
#include "../inc/exch.h"
void sim::stepchange()
{
    config::printline(config::Info);

    libconfig::Setting &setting = config::cfg.lookup("stepchange");
    double Tstart=0.0,Tfinal=0.0,met=0.0,et=0.0,rate=0.0;
    setting.lookupValue("T_start",Tstart);
    FIXOUT(config::Info,"Initial temperature:" << Tstart << std::endl);
    setting.lookupValue("T_final",Tfinal);
    FIXOUT(config::Info,"Final temperature:" << Tfinal << std::endl);
    setting.lookupValue("Rate",rate);
    FIXOUT(config::Info,"Rate of cooling/heating:" << rate << " K/s" << std::endl);
    setting.lookupValue("RunTime",met);
    setting.lookupValue("EquilTime",et);
    FIXOUT(config::Info,"Equilibration time:" << et << " seconds" << std::endl);
    FIXOUT(config::Info,"Run time:" << et << " seconds" << std::endl);

    unsigned int mrts=int(met/llg::dt),ets=int(et/llg::dt);
    std::string opf;
    bool outmag=false;
    std::ofstream magout,emagout;
    std::string MagFilestr;
    setting.lookupValue("outmag",outmag);
    FIXOUT(config::Info,"Output magnetization data for each temperature:" << config::isTF(outmag) << std::endl);
    if(outmag==true)
    {
        setting.lookupValue("MagFileName",MagFilestr);
        FIXOUT(config::Info,"Magnetization files output to:" << MagFilestr << std::endl);
    }

    config::printline(config::Info);
    llg::T=Tstart;
    if(outmag)
    {
        std::stringstream MagFilesstr,eqsstr;
        MagFilesstr << MagFilestr << "_" << Tstart << "_" << Tfinal << ".dat";
        std::string mopf=MagFilesstr.str();
        magout.open(mopf.c_str());
        if(!magout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for outputting magnetization data");
        }
    }

    for(unsigned int t = 0 ; t < ets ; t++)
    {
        if(t%spins::update==0)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);
            magout << Tstart << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            if(t%(spins::update*100)==0)
            {
//                util::outputSpinsVTU(t);
            }
        }

        llg::integrate(t);
    }
    double Ti=Tstart;
    double Tip1=0.0;
    for(unsigned int t = ets ; t < mrts+ets ; t++)
    {
        Tip1=rate*llg::dt+Ti;
        llg::T=Tip1;
        Ti=Tip1;
        if(rate<0 && Ti<Tfinal)
        {
            Ti=Tfinal;
            Tip1=Tfinal;
            llg::T=Tfinal;
        }
        if(rate > 0 && Ti>Tfinal)
        {
            Ti=Tfinal;
            Tip1=Tfinal;
            llg::T=Tfinal;
        }
		if(t%spins::update==0)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);

            magout << llg::T << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
        }
        llg::integrate(t);

    }

    magout.close();
    if(magout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close file for outputting magnetization data.");
    }
}
