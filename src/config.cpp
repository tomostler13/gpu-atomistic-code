// File: config.cpp
// Author:Tom Ostler
// Last-modified: 20 Feb 2013 12:47:27
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdarg>
#include "../inc/error.h"
#include "../inc/random.h"
#include "../inc/util.h"
#include <cassert>
#include "../inc/config.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace config
{
    libconfig::Config cfg;
    bool lcf=false;
    bool incdip=true;
    unsigned int seed=0;
    std::ofstream Info;
    void initConfig(int argc,char *argv[])
    {
        libconfig::Config cfg;
        if(argc < 2)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You must give a config file, exiting");
        }
        try
        {
            cfg.readFile(argv[1]);
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
        //for getting the date and time
        time_t now = time(0);
        char *dtime=ctime(&now);
        std::string iffstr=cfg.lookup("OutputFile");
        seed = cfg.lookup("seed");
        incdip=cfg.lookup("include_dipolar");
        FIXOUT(config::Info,"Including dipolar terms:" << isTF(incdip) << std::endl);
        Random::seed(seed,seed+100);
        Info.open(iffstr.c_str());
        //open the output info file
        if(!Info.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open info file for writing to, exiting");
        }
        else
        {
			Info.width(45);Info << std::right << "*" << "**Build details***" << std::endl;
            Info.width(75);Info << std::left << "Opened output file (this file): " << iffstr << std::endl;
        }
        Info.width(75);Info << std::left << "Run Date/time:" << dtime << std::flush;
		if(GITDIRTY!="0")
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errWarning("Warning: Your git build is dirty, you should not use this code for production");
			FIXOUT(Info,"Git SHA:" << GIT_SHA1 << ", Dirty" << std::endl);
		}
		else
		{
			FIXOUT(Info,"Git SHA:" << GIT_SHA1 << ", Clean" << std::endl);
		}
		FIXOUT(Info,"Compile Data/Time:" << __DATE__ << ", " << __TIME__ << std::endl);
		FIXOUT(Info,"Compiler:" << COMP << std::endl);
		FIXOUT(Info,"Compiled on machine:" << HOSTNAME << std::endl);
        FIXOUT(Info,"Localhost:" << util::exec("hostname") << std::endl);
        FIXOUT(Info,"Seed:" << seed << std::endl);

        assert(seed>0);
        lcf=true;
    }
    void printline(std::ofstream& os)
    {
        for(unsigned int i = 0 ; i < 125 ; i++)
        {
            os << "-";
        }
        os << std::endl;
    }
    std::string isTF(bool cb)
    {
        if(cb)
        {
            return("TRUE");
        }
        else
        {
            return("FALSE");
        }
    }

}
