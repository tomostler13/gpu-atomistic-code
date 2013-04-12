// File: test_nanowire.cpp
// Author: Tom Ostler
// Created: 10 April 2013
// Last-modified: 12 Apr 2013 18:58:30
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
void sim::test_nanowire(int argc,char *argv[])
{
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Test Nanowire details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("test_nanowire");

    //low temperature reservoir
    double LTR=0.0;
    //high temperature reservoir
    double HTR=0.0;
    //temperature of wire
    double TOW=0.0;
    //equilibration and runtime
    double et=0.0,rt=0.0;
    //steps
    bool outmag=false;
    std::string MagFilestr;
    setting.lookupValue("EquilTime",et);
    setting.lookupValue("RunTime",rt);
    setting.lookupValue("Wire_Temperature",TOW);
    setting.lookupValue("High_Temperature",HTR);
    setting.lookupValue("Low_Temperature",LTR);
    setting.lookupValue("outmag",outmag);

    unsigned int ets=int(et/llg::dt),rts=int(rt/llg::dt);
    FIXOUT(config::Info,"Output magnetization data for each temperature:" << config::isTF(outmag) << std::endl);
    if(outmag==true)
    {
        setting.lookupValue("MagFileName",MagFilestr);
        FIXOUT(config::Info,"Magnetization files output to:" << MagFilestr << std::endl);
    }
    std::ofstream magout;
    if(outmag)
    {
        std::stringstream MagFilesstr,eqsstr;
        MagFilesstr << MagFilestr << "_" << LTR << "_" << TOW << "_" << HTR << ".dat";
        std::string mopf=MagFilesstr.str();
        magout.open(mopf.c_str());
        if(!magout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for outputting magnetization data");
        }
    }

    FIXOUT(config::Info,"Temperature of high temperature reservoir:" << HTR << std::endl);
    FIXOUT(config::Info,"Temperature of low temperature reservoir:" << LTR << std::endl);
    FIXOUT(config::Info,"Temperature of wire:" << TOW << std::endl);
    FIXOUT(config::Info,"Equilibration time:" << et << std::endl);
    FIXOUT(config::Info,"Run time:" << rt << std::endl);
    FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    FIXOUT(config::Info,"Number of runtime steps:" << rts << std::endl);
    FIXOUT(config::Info,"Setting temperature:" << std::flush);
    std::cout << llg::osT.size() << std::endl;
    for(unsigned int i = 0 ; i < geom::nspins ; i++)
    {
        int coords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
        std::cout << coords[2] << std::endl;
        if(coords[2]<geom::cut0)
        {
            llg::osT[i]=LTR;
            std::cout << "BALSDJKLASDFJKLSADF" << std::endl;
        }
        else if(coords[2]>=geom::cut0 && coords[2]<geom::cut1)// < geom::cut1 && coords[0] > (geom::dim[0]-geom::width)/2 && coords[0] < (geom::dim[0]+geom::width)/2 && coords[1] > (geom::dim[1]-geom::width)/2 && coords[1] < (geom::dim[1]+geom::width)/2)
        {
            llg::osT[i]=TOW;
            std::cout << "________________________________" << std::endl;
        }
        else if(coords[2]>=geom::cut1)
        {
            llg::osT[i]=HTR;
            std::cout << "Setting" << std::endl;
        }
    }
    std::cin.get();
    SUCCESS(config::Info);
    if(llg::rscf)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("You have set output of correlation function to true when you do not have translational invariance, setting to false.");
        llg::rscf=false;
    }
    std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
    for(unsigned int t = 0 ; t < ets ; t++)
    {
        llg::integrate(t,llg::osT,exch::xadj,exch::adjncy);
        if(t%spins::update==0)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);
            magout << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            util::outputSpinsVTU(t);
        }
    }

    for(unsigned int t = ets ; t < rts+ets ; t++)
    {
        llg::integrate(t,llg::osT,exch::xadj,exch::adjncy);
        if(t%spins::update==0)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);
            magout << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            util::outputSpinsVTU(t);
        }
    }



}
