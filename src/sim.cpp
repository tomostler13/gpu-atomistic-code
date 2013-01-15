// File: sim.cpp
// Author:Tom Ostler
// Last-modified: 11 Dec 2012 15:29:23
#include "../inc/config.h"
#include "../inc/sim.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/tdp.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>
namespace sim
{
    bool simi=false;
    double dt=0;
    double rdt=0;
    double llbpf=0;
    std::string sim_type;
    void initSim(int argc,char *argv[])
    {
        assert(mat::mi);
        assert(geom::gi);
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
            setting.lookupValue("dt",dt);
            setting.lookupValue("sim_type",sim_type);
        }
        catch(const libconfig::SettingTypeException &stex)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting type error");
        }
        FIXOUT(config::Info,"Simulation selected:" << sim_type << std::endl);
        FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
        rdt=dt*mat::gamma;
        //sigma does not change it is equal to sqrt(2*|gamma|*k_B/(dt*Ms*V))
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            tdp::sysSigma[i]=sqrt(1.38e-23*2.0/(rdt*mat::Ms*geom::gsV));
        }

        llbpf=-1;
        simi=true;
    }

}
