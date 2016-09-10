// File: sim.cpp
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 10 Sep 2016 18:47:39
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
        if(!config::cfg.exists("sim")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting sim does not exist. It is required.");
        }
        libconfig::Setting &setting = config::cfg.lookup("sim");
        if(!setting.lookupValue("sim_type",sim_type)0
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Error reading simulation type (sim.sim_type (string)). Options are mvt, laserheating");
        }
        FIXOUT(config::Info,"Simulation selected:" << sim_type << std::endl);
    }

}
