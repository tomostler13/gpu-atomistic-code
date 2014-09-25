// File: sim.cpp
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 24 Sep 2014 15:03:20
#include "../inc/config.h"
#include "../inc/sim.h"
#include "../inc/geom.h"
#include "../inc/defines.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>
namespace sim
{
    std::string sim_type;
    void initSim(int argc,char *argv[])
    {
        if(argc < 2)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You must give a config file, exiting");
        }
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

        libconfig::Setting &setting = config::cfg.lookup("sim");
        try
        {
            setting.lookupValue("sim_type",sim_type);
        }
        catch(const libconfig::SettingTypeException &stex)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting type error");
        }
        FIXOUT(config::Info,"Simulation selected:" << sim_type << std::endl);
    }

}
