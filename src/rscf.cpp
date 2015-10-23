// File: sf_glob.cpp
// Note: originall dsf_glob.cpp
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 23 Oct 2015 17:33:45
#include "../inc/llg.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include "../inc/spins.h"
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

//Reads the parameters and sets up parts of the SF calculation
namespace rscf
{
    //calculate correlation function?
    bool ccf=false;
    //number of points
    unsigned int nr=0;
    Array2D<bool> cab;
    Array2D<fftw_plan> rscfP;
    Array2D<unsigned int> rpoints;
    void initRSCF(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Correlation function details***" << std::endl;
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
        bool errstatus=false;
        std::string cffile;
        libconfig::Setting &setting = config::cfg.lookup("rscf");
        errstatus=setting.lookupValue("Calculate",ccf);
        //do some error handling
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read whether the correlation function calculation is selected.");
        }
        else
        {
            FIXOUT(config::Info,"Outputting correlation function?:" << config::isTF(ccf) << std::endl);
            if(ccf==false)
            {
                //we do not need to any more processing of the dsf information
                return;
            }
            errstatus=setting.lookupValue("InputFile",cffile);
            if(errstatus==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You specified the use of the real space correlation function calculation, however you have to speficy an input file (rscf:InputFile) in the config file.");
            }
            else
            {
                FIXOUT(config::Info,"Correlation function input file:" << cffile << std::endl);
            }
        }
        if(ccf)
        {
            try
            {
                config::cfg.readFile(cffile.c_str());
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
            libconfig::Setting &nsetting = config::cfg.lookup("rscf");
            cab.resize(3,3);
            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                try
                {
                    cab(0,i)=nsetting["Cx_beta"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Setting " << snf.getPath();
                    error::errMessage(errsstr);
                }
                try
                {
                    cab(1,i)=nsetting["Cy_beta"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsste << "Setting not found exception caught. Setting " << snf.getPath();
                    error::errMessage(errsstr);
                }
                try
                {
                    cab(2,i)=nsetting["Cz_beta"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsste << "Setting not found exception caught. Setting " << snf.getPath();
                    error::errMessage(errsstr);
                }
            }

        }
    }
}
