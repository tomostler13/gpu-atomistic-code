// File: fmr.cpp
// Author: Tom Ostler
// Created: 14 June 2013
// Last-modified: 14 Jun 2013 09:37:23
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
void sim::fmr(int argc,char *argv[])
{
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**FMR details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("fmr");
    double redfreq=0.0,Bdrive=0.0,Ncycles=0,MinCycles=0,MaxCycles=0;
    double temp=0.0,ConvVar=2e-11,ConvMean=1e-8;

    std::string opfs;
    try
    {
        setting.lookupValue("temperature",temp);
        setting.lookupValue("redfreq",redfreq);
        setting.lookupValue("Bdrive",Bdrive);
        setting.lookupValue("ConvVar",ConvVar);
        setting.lookupValue("ConvMean",ConvMean);
        setting.lookupValue("Ncycles",Ncycles);
        setting.lookupValue("MaxCycles",MaxCycles);
        setting.lookupValue("MinCycles",MinCycles);
        setting.lookupValue("MagFile",opfs);
    }
    catch(const libconfig::SettingTypeException &stex)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting type error");
    }
    if(redfreq < 1e-12)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Frequency is too small, make sure it is read in properly");
    }

    FIXOUT(config::Info,"Temperature:" << temp << std::endl);
    FIXOUT(config::Info,"Driving frequency/gamma:" << redfreq << std::endl);
    FIXOUT(config::Info,"Magnitude of driving field:" << Bdrive << " [T]" << std::endl);
    FIXOUT(config::Info,"Number of equilibration cycles:" << Ncycles << std::endl);
    FIXOUT(config::Info,"Minimum cycles to average over:" << MinCycles << std::endl);
    FIXOUT(config::Info,"Maximum cycles to average over:" << MaxCycles << std::endl);
    FIXOUT(config::Info,"Average tolerance:" << ConvMean << std::endl);
    FIXOUT(config::Info,"Variance tolerance:" << ConvVar << std::endl);
    FIXOUT(config::Info,"Outputting magnetization data to:" << opfs << std::endl);
    std::ofstream emagfile(opfs.c_str());
    if(!emagfile.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file for writing magnetization data");
    }
    unsigned int ntspc=0;//number of timesteps per cycles
    //----------------------find adjusted frequency----------------------------------------------
    {//only want these variables to have temporary scope
        config::Info << std::setprecision(10);
        double Tp=2.0*M_PI/(redfreq*mat::gamma);//time period
        double nts=Tp/llg::dt;//the correct number of timesteps may not be
        FIXOUT(config::Info,"Ideal number of timesteps per cycle:" << nts << std::endl);
        //an integer
        ntspc=int(nts+0.5);//round to nearest integer
        FIXOUT(config::Info,"Actual number of timesteps per cycle:" << ntspc << std::endl);
        double newfreq=2.0*M_PI/(double(ntspc)*llg::dt);
        FIXOUT(config::Info,"New frequency with corrected timesteps:" << newfreq << " rad/s" << std::endl);
        newfreq=newfreq/mat::gamma;
        redfreq=newfreq;
        FIXOUT(config::Info,"New reduced frequency:" << redfreq << std::endl);
    }
}
