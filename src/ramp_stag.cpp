// File: ramp_field.cpp
// Author: Tom Ostler
// Created: 15 May 2023
// Last-modified: 31 May 2023 17:03:08
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
void sim::ramp_stag()
{
    //rampinc is the increment (from 0 to 1 of the original staggered field value)
    double ramptime=0.0,runtime=0.0;
    // How much the staggered field will change at each timestep
    double dHstg=0.0;
    //rf is a unit vector to specify the direction in which the field
    //is changed
    unsigned int rampts=0,runts=0;
    std::string filestr;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Staggered Field Ramp Details***" << std::endl;
    //The k-vector information is read in the file sf.cpp and sf_glob.cpp

    if(!config::cfg.exists("ramp_stag"))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting \"ramp_stag\" does not exist. Check your config file.");
    }
    libconfig::Setting &setting = config::cfg.lookup("ramp_stag");
    bool errstatus=false;
    if(setting.lookupValue("RampTime",ramptime))
    {
        FIXOUT(config::Info,"Ramping time:" << ramptime << " [s]" << std::endl);
        rampts = static_cast<unsigned int>(ramptime/llg::dt+0.5);
        FIXOUT(config::Info,"Number of ramping timesteps:" << rampts << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting ramp_stag.RampTime (double)");
    }

    rampinc=1.0/static_cast<double>(rampts);
    FIXOUT(config::Info,"Increment of ramp (from 0 to 1):" << rampinc << std::endl);

    if(setting.lookupValue("RunTime",runtime))
    {
        FIXOUT(config::Info,"Running time:" << runtime << " [s]" << std::endl);
        runts = static_cast<unsigned int>(runtime/llg::dt+0.5);
        FIXOUT(config::Info,"Number of running timesteps:" << runts << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting ramp_stag.RunTime (double)");
    }

    errstatus=setting.lookupValue("Temperature",llg::T);
    if(errstatus)
    {
        FIXOUT(config::Info,"Temperature:" << llg::T << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read temperature for ramp_stag (ramp_stag:Temperature (double))");
    }

    double rv=0.0;
    for(unsigned int t = 0 ; t < rampts ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);

        }
        llg::integrate(t);
        fields::incStagFields(rv);
        rv+=rampinc;
    }

    for(unsigned int t = rampts ; t < rampts+runts ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);

        }
        llg::integrate(t);
    }

}
