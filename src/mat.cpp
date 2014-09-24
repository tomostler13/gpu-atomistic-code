// File: mat.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 24 Sep 2014 10:31:56
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
   }
}
