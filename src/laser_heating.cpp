// File: laser_heating.cpp
// Author: Tom Ostler
// Created: 24 Nov 2014
// Last-modified: 24 Nov 2014 13:13:51
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
//function prototype
double pump_power(double);
void sim::laser_heating(int argc,char *argv[])
{
    bool opsf=false;//OutPut Structure Factor info
    unsigned int ets = 0 , rts = 0 , num_pulses=0;
    double et = 0.0 , rt = 0.0 , gamma_e = 0.0 , Cl = 0.0 , G_ep = 0.0 , Pump_Time = 50e-15, ct = 0.0 , T_el_max = 0;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Laser Heating Simulation details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("laserheating");
    bool errstatus=setting.lookupValue("EquilibrationTime",et);
    if(errstatus)
    {
        FIXOUT(config::Info,"Equilibration time:" << et << " [s]" << std::endl);
        ets=static_cast<unsigned int>((et/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    }
    errstatus=setting.lookupValue("RunTime",rt);
    if(errstatus)
    {
        FIXOUT(config::Info,"Run time:" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run:" << rts << std::endl);
    }
    errstatus=setting.lookupValue("InitialTemperature",llg::T);
    if(errstatus)
    {
        FIXOUT(config::Info,"Initial Temperature:" << llg::T << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read temperature (laserheating:InitialTemperature (double))");
    }
    errstatus=setting.lookupValue("Cl",Cl);
    if(errstatus)
    {
        FIXOUT(config::Info,"Lattice specific heat:" << Cl << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read lattice specific heat (laserheating:Cl (double))");
    }
    errstatus = setting.lookupValue("PumpTime",Pump_Time);
    if(errstatus)
    {
        FIXOUT(config::Info,"Laser pump time (width):" << Pump_Time << " [s]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read laser pump time (laserheating:PumpTime (double)).");
    }
    errstatus = setting.lookupValue("Gamma_e",gamma_e);
    if(errstatus)
    {
        FIXOUT(config::Info,"gamma_e (Ce(Te)=gamma_e*Te):" << gamma_e << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the gamma_e (Ce(T)=gamma_e*Te) (laserheating:Gamma_e (double))");
    }
    errstatus = setting.lookupValue("CoolingTime",ct);
    if(errstatus)
    {
        FIXOUT(config::Info,"Cooling time:" << ct << " [s]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the cooling time (laserheating:CoolingTime (double))");
    }
    errstatus = setting.lookupValue("OutputStructureFactor",opsf);
    if(errstatus)
    {
        FIXOUT(config::Info,"Output the information for the structure factor?" << config::isTF(opsf) << std::endl);
        if(opsf)
        {
            if(sf::csf==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("The flag laserheating:OutputStructureFactor was set to true and you have not set the structure factor and k-vectors up (sf:Calculate)");
            }
        }
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the state of the outputting of the structure factor (laserheating:OutputStructureFactor (bool))");
    }
    errstatus = setting.lookupValue("MaxElectronTemperature",T_el_max);
    if(errstatus)
    {
        FIXOUT(config::Info,"Input max electron temperature:" << T_el_max << std::endl);
        //convert to power
        T_el_max = pump_power(T_el_max);
        FIXOUT(config::Info,"Laser power:" << T_el_max << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the max electron temperature (laserheating:MaxElectronTemperature (double))");
    }
    errstatus = setting.loopkupValue("NumPulses",num_pulses);
    if(errstatus)
    {
        FIXOUT(config::Info,"Number of pulses:" << num_pulses << std::endl);
    }
    error::errPreamble(__FILE__,__LINE__);
    double pulse_delays[num_pulses];
    FIXOUT(config::Info,"Pulse delays:"<< " [" << std::flush);
    for(unsigned int i = 0 ; i < num_pulses ; i++)
    {
        pulse_delays[i]=setting["PulseDelays"][i];
        if(i<(num_pulses-1))
        {
            config::Info << pulse_delays[i] << " , ";
        }
    }
    config::Info << pulse_delays[num_pulses-1] << " ]" << std::endl;
    double pulse_scale[num_pulses];
    FIXOUT(config::Info,"Pulse scaling factors:"<< " [" << std::flush);
    for(unsigned int i = 0 ; i < num_pulses ; i++)
    {
        pulse_scale[i]=setting["PulseScale"][i];
        if(i<(num_pulses-1))
        {
            config::Info << pulse_scale[i] << " , ";
        }
    }
    config::Info << pulse_scale[num_pulses-1] << " ]" << std::endl;

}

double pump_power(double x)//we do not pass by reference here deliberately because T_el_max set
{
    return((-5.9408+0.00722*x+4.11859e-5*(x*x)-3.66025e-10*x*x*x)*1.152E20);
}
