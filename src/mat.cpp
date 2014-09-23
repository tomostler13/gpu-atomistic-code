// File: mat.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 23 Sep 2014 13:13:14
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
        //resize all of the arrays that are number of magetic species (geom::nms) long
        lambda.resize(geom::nms);gamma.resize(geom::nms);mu.resize(geom::nms);sigma.resize(geom::nms);
        for(unsigned int i = 0 ; i < geom::nms ; i++)
        {
            //get damping constants
            lambda[i]=setting["lambda"][i];
            std::stringstream sstr;
            sstr << "Coupling constants (lambda) for species " << i << ":";
            std::string str=sstr.str();
            FIXOUT(config::Info,str.c_str() << lambda[i] << std::endl);
            //get gyromagnetic ratios
            gamma[i]=setting["gamma"][i];
            sstr.str("");
            sstr << "Gyromagnetic ratio for species " << i << ":";
            str=sstr.str();
            FIXOUT(config::Info,str.c_str() << gamma[i] << " [in units of gyromagnetic ratio]\n");
            //get the magnetic moments
            mu[i]=setting["mu"][i];
            sstr.str("");sstr << "Magnetic moment for species " << i << ":";
            str=sstr.str();
            FIXOUT(config::Info,str.c_str() << mu[i] << " [muB]" << std::endl);
        }
    }
}
