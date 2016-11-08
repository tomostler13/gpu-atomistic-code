// File: sim.cpp
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 18 Oct 2016 13:20:00
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
    void initSim()
    {
        if(!config::cfg.exists("sim"))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting sim does not exist. It is required.");
        }
        libconfig::Setting &setting = config::cfg.lookup("sim");
        if(!setting.lookupValue("sim_type",sim_type))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Error reading simulation type (sim.sim_type (string)). Options are mvt, laserheating");
        }
        FIXOUT(config::Info,"Simulation selected:" << sim_type << std::endl);
    }

}
