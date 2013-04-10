// File: stepchange.cpp
// Author: Tom Ostler
// Created: 29 Mar 2013
// Last-modified: 10 Apr 2013 13:52:50
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
void sim::stepchange(int argc,char *argv[])
{
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Step Change details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("stepchange");
    double Tstart=0.0,Tfinal=0.0,met=0.0,et=0.0;
    setting.lookupValue("T_start",Tstart);
    FIXOUT(config::Info,"Initial temperature:" << Tstart << std::endl);
    setting.lookupValue("T_final",Tfinal);
    FIXOUT(config::Info,"Final temperature:" << Tfinal << std::endl);
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
        llg::integrate(t);
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
    }
    llg::T=Tfinal;
    for(unsigned int t = ets ; t < mrts+ets ; t++)
    {
        llg::integrate(t);
        if(t%spins::update==0)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);

            magout << Tfinal << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            if(t%(spins::update*100)==0)
            {
//                util::outputSpinsVTU(t);
            }
        }
    }

    magout.close();
    if(magout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close file for outputting magnetization data.");
    }
}
