// File: thermal_hyst.cpp
// Author: Tom Ostler
// Created: 7th June 2015
// Last-modified: 07 Jun 2015 16:21:49
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>
#include "../inc/error.h"
#include "../inc/arrays.h"
#include "../inc/defines.h"
#include "../inc/util.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/intmat.h"
#include "../inc/fields.h"
#include "../inc/llg.h"
#include "../inc/sim.h"
#include "../inc/sf.h"
void sim::thermal_hyst(int argc,char *argv[])
{
    double Tstart=0.0;
    double dT=0.0;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Thermal Hysteresis details***" << std::endl;
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
    bool checkset=config::cfg.exists("thermal_hyst");
    if(!checkset)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting \"thermal_hyst\" not found. Check your config file.");
    }
    libconfig::Setting &setting=config::cfg.lookup("thermal_hyst");
    if(setting.lookupValue("StartTemp",Tstart))
    {
        FIXOUT(config::Info,"Starting temperature:" << Tstart << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting thermal_hyst.StartTemp (double)");
    }
    if(setting.lookupValue("IncTemp",dT))
    {
        FIXOUT(config::Info,"Temperature increment:" << dT << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting thermal_hyst.IncTemp (double)");
    }
    if(setting.lookupValue("EquilTime",et))
    {
        FIXOUT(config::Info,"Equilibration time per temperature:" << et << " [s]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting thermal_hyst.EquilTime (double)");
    }
    if(setting.lookupValue("RunTime",rt))
    {
        FIXOUT(config::Info,"Run time per temperature:" << rt << " [s]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Coudl not read the setting thermal_hyst.RunTime (double)");
    }
}
