// File: mat.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 21 Jan 2013 14:45:11
#include "../inc/mat.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/arrays.h"
#include <libconfig.h++>
#include <cassert>
#include <sstream>
#include <string>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace mat
{
    //bath coupling (no units)
    double lambda=0.0;
    //Gyromagnetic ratio
    double gamma=0.0;
    //Bohr magneton
    double muB=9.27e-24;
    //magnetic moment (muB)
    double mu=2.0;
    //thermal term prefactor
    double sigma = 0.0;
    void initMat(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Material details***" << std::endl;
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

        libconfig::Setting &setting = config::cfg.lookup("mat");

        setting.lookupValue("lambda",lambda);
        FIXOUT(config::Info,"Coupling constant (lambda):" << lambda << std::endl);
        setting.lookupValue("gamma",gamma);
        FIXOUT(config::Info,"Gyromagnetic ratio:" << gamma << " T^-1 s^-1\n");
        muB=9.27e-24;
        setting.lookupValue("mu",mu);
        FIXOUT(config::Info,"Magnetic moment:" << mu << " muB" << std::endl);
    }
}
