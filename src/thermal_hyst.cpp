// File: thermal_hyst.cpp
// Author: Tom Ostler
// Created: 7th June 2015
// Last-modified: 27 Nov 2015 19:01:51
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
    bool outputSpinConfig=false;
    bool checkset=config::cfg.exists("thermal_hyst");
    if(!checkset)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting \"thermal_hyst\" not found. Check your config file.");
    }
    unsigned int numint=0,ets=0;
    double et=0;
    unsigned int sof=0;
    libconfig::Setting &setting=config::cfg.lookup("thermal_hyst");
    if(!setting.lookupValue("OutputSpinConfig",outputSpinConfig))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting thermal_hyst.OutputSpinConfig (bool)");
    }

    FIXOUT(config::Info,"Output spin config?:" << config::isTF(outputSpinConfig) << std::endl);
    if(outputSpinConfig)
    {
        if(!setting.lookupValue("SpinOutputFreq",sof))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You selected to output the spin map but the code could not read the setting thermal_hyst.SpinOutputFreq (unsigned int)");
        }
        else
        {
            FIXOUT(config::Info,"Spin update (in units of spin update):" << sof << " [timesteps]" << std::endl);
        }
    }
    if(setting.lookupValue("EquilTime",et))
    {
        FIXOUT(config::Info,"Equilibration time:" << et << std::endl);
        ets = static_cast<unsigned int>(et/llg::dt+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting thermal_hyst.EquilTime (double)");
    }
    if(setting.lookupValue("NumRateIntervals",numint))
    {
        FIXOUT(config::Info,"Number of rate intervals:" << numint << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting thermal_hyst.NumRateIntervals (int)");
    }
    Array<double> temps,rates;
    temps.resize(numint+1);
    rates.resize(numint);
    temps.IFill(0);
    rates.IFill(0);
    for(unsigned int i = 0 ; i < numint+1 ; i++)
    {
        try
        {
            temps[i]=setting["Temps"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
    }
    Array<double> dT;
    dT.resize(numint);
    dT.IFill(0);
    for(unsigned int i = 0 ; i < numint ; i++)
    {
        std::stringstream sstr,sstr1;
        sstr << "Interval " << i << " has temperature range ";
        sstr1 << temps[i] << " -> " << temps[i+1] << " [K]";
        std::string str=sstr.str(),str1=sstr1.str();
        FIXOUT(config::Info,str << str1 << std::endl);
        try
        {
            rates[i]=setting["Rates"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
        std::stringstream sstrn;
        sstr << "Interval " << i << " has rate:";
        std::string strn=sstrn.str();
        FIXOUT(config::Info,strn << rates[i] << std::endl);
        dT[i]=rates[i]*llg::dt;
        FIXOUT(config::Info,"Temperature increase (per update):" << dT[i] << " [K]" << std::endl);

    }
    std::ofstream Ttime("T_time.dat");
    llg::T=temps[0];
    if(!Ttime.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file T_time.dat");
    }
    Ttime << std::setprecision(6);
    for(unsigned int t = 0 ; t < ets ; t++)
    {
        if(t%spins::update==0)
        {
            Ttime << static_cast<double>(t)*llg::dt << "\t" << llg::T << std::endl;
            util::calc_mag();
            util::output_mag(t);
        }
        llg::integrate(t);
    }
    unsigned int time_counter=ets;
    if(rates[0]>0)//then we are heating
    {
        for(unsigned int i = 0 ; i < numint ; i++)
        {
            for(llg::T = temps[i] ; llg::T < temps[i+1] ; llg::T+=dT[i])
            {
                if(time_counter%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(time_counter);
                    Ttime << static_cast<double>(time_counter)*llg::dt << "\t" << llg::T << std::endl;
                }

                llg::integrate(time_counter);
                time_counter++;
            }
        }
    }
    else//we are cooling
    {
        for(unsigned int i = 0 ; i < numint ; i++)
        {
            for(llg::T = temps[i] ; llg::T > temps[i+1] ; llg::T+=dT[i])
            {
                if(time_counter%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(time_counter);
                    Ttime << static_cast<double>(time_counter)*llg::dt << "\t" << llg::T << std::endl;
                }

                llg::integrate(time_counter);
                time_counter++;
            }
        }
    }
    Ttime.close();
    if(Ttime.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close file T_time.dat");
    }
}
