// File: config.cpp
// Author:Tom Ostler
// Last-modified: 07 Oct 2014 10:43:11
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
#include "../inc/config.h"
#include "../inc/defines.h"
#include <cassert>
namespace config
{
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
        char c[] = "hostname";
        FIXOUT(Info,"Localhost:" << util::exec(c) << std::endl);
        FIXOUT(Info,"Seed:" << seed << std::endl);

        libconfig::Setting &setting = cfg.lookup("system");
        setting.lookupValue("include_dipole",inc_dip);
        setting.lookupValue("Exchange_method",exchmeth);
        std::cout << "Excnage  method read as\t" << exchmeth << std::endl;
        setting.lookupValue("Dipole_method",dipmeth);
        //read whether we want PBC's or not
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            pbc[i]=setting["Periodic_Boundaries"][i];
        }
        FIXOUTVEC(config::Info,"Periodic boundary conditions",isTF(pbc[0]),isTF(pbc[1]),isTF(pbc[2]));
        //So that we are not comparing strings all over the place
        //we assign a method to an integer. This also ensures we have
        //written the correct string
        if(exchmeth=="fft")
        {
            exchm=0;
        }
        else if(exchmeth=="DIA")
        {
            exchm=1;
        }
        else if(exchmeth=="CSR")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CSR Coming soon");
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Exchange interaction method (Exchange_method) not recognized, please select either fft, DIA or CSR");
        }
        if(dipmeth=="fft")
        {
            dipm=0;
        }
        else if(dipmeth=="DIA")
        {
            dipm=1;
        }
        else if(dipmeth=="CSR")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("CSR Coming soon");
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Dipole interaction method (Dipole_method) not recognized, please select either fft, DIA or CSR");
        }
        //check if you are using a sparse matrix multiplication (somewhere)
        //if so we have the option of not calculating the off-diagonals
        if((dipm+exchm)>0)
        {
            setting.lookupValue("Off_Diagonals",offdiag);
            FIXOUT(config::Info,"Include off-diagonals in sparse matrix multiplication:" << config::isTF(offdiag) << std::endl);
        }
        //if we are using the fft at the moment we cannot have PBC's
        if((dipm==0 || exchm==0) && (pbc[0]==true || pbc[1]==true || pbc[2]))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You cannot currently use the fft method for calculating exchange or dipole-dipole fields and have periodic boundary conditions.\nIf you want to use PBC's then select a matrix multiplication method.");
        }
        FIXOUT(config::Info,"Exchange field calculation method:" << exchmeth << std::endl);
        FIXOUT(config::Info,"Dipole-dipole field calculation method:" << dipmeth << std::endl);
        assert(seed>0);
        lcf=true;
    }
    void openLogFile()
    {
        if(!Log.is_open())
        {
            Log.open("log.dat");
            if(!Log.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open log file");
            }
        }
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
