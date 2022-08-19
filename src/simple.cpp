// File: simple.cpp
// Author: Tom Ostler
// Created: 19 August 2022
// Last-modified: 19 Aug 2022 13:36:28
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
void sim::simple()
{
    unsigned int num_samples=0;
    double et=0.0,rt=0.0;
    //ets = equilibration time
    //rts = run timesteps
    unsigned int ets=0,rts=0;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Simple calculation details***" << std::endl;
    //The k-vector information is read in the file sf.cpp and sf_glob.cpp
    if(!config::cfg.exists("simple"))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting \"simple\" does not exist. Check your config file.");
    }
    if(!config::cfg.exists("simple"))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Simple settings cannot be found, check your config file.");
    }
    libconfig::Setting &setting = config::cfg.lookup("simple");
    if(setting.lookupValue("EquilibrationTime",et))
    {
        FIXOUT(config::Info,"Equilibration time:" << et << " [s]" << std::endl);
        ets=static_cast<unsigned int>((et/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not find the equilibration time, defaulting to 15ps");
        et=15e-12;
        FIXOUT(config::Info,"Equilibration time (default):" << et << " [s]" << std::endl);
        ets=static_cast<unsigned int>((et/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps (default):" << ets << std::endl);
    }
    if(setting.lookupValue("RunTime",rt))
    {
        FIXOUT(config::Info,"Run time:" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run:" << rts << std::endl);
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not find the run time, defaulting to 100ps");
        rt=100e-12;
        FIXOUT(config::Info,"Run time (default):" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run (default):" << rts << std::endl);
    }
    if(setting.lookupValue("Temperature",llg::T))
    {
        FIXOUT(config::Info,"Temperature:" << llg::T << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read temperature for simple (simple:Temperature (double)). Must be set.");
    }


    fields::setStagZero();

    for(unsigned int t = 0 ; t < ets ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
        }
        llg::integrate(t);
    }

    fields::setStagField();

    for(unsigned int t = ets ; t < (rts+ets) ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
        }

        llg::integrate(t);
    }
    unsigned int tend=-1;
    util::outputSpinsVTU(tend);


    //OUTPUT DW PROFILE AT END
    std::ofstream dwofs("dw_grid.out");
    if(dwofs.is_open())
    {
        dwofs << geom::nspins << std::endl;
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            dwofs << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << spins::Sx[i] << "\t" << spins::Sy[i] << "\t" << spins::Sz[i] << std::endl;
        }
        dwofs.close();
        if(dwofs.is_open())
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Couldn't close dw_grid.out file");
        }
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Couldn't open dw_grid.out file, redirection to std::cout, this could get messy.");
        std::cout << geom::nspins << std::endl;
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << spins::Sx[i] << "\t" << spins::Sy[i] << "\t" << spins::Sz[i] << std::endl;
        }

    }

}
