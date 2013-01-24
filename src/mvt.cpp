// File: mvt.h
// Author: Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 23 Jan 2013 22:28:52
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
void sim::MvT(int argc,char *argv[])
{
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**M(T) details***" << std::endl;
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

        libconfig::Setting &setting = config::cfg.lookup("mvt");
		double lT=0.0,uT=0.0,dT=0.0,convmean=0.0,convvar=0.0;
		setting.lookupValue("lower_temp",lT);
		FIXOUT(config::Info,"Lower temperature:" << lT << std::endl);
		setting.lookupValue("upper_temp",uT);
		FIXOUT(config::Info,"Upper temperature:" << uT << std::endl);
		setting.lookupValue("temp_step",dT);
		FIXOUT(config::Info,"Temperature step:" << dT << std::endl);
		setting.lookupValue("mean_tolerance",convmean);
		setting.lookupValue("variance_tolerance",convvar);
		FIXOUT(config::Info,"Converging mean to:" << convmean << std::endl);
		FIXOUT(config::Info,"Converging variance to:" << convvar << std::endl);
		util::RunningStat MS;


}
