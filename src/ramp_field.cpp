// File: ramp_field.cpp
// Author: Tom Ostler
// Created: 13 May 2015
// Last-modified: 12 Aug 2015 11:27:51
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>
#include "../inc/error.h"
#include "../inc/arrays.h"
#include "../inc/defines.h"
#include "../inc/util.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/intmat.h"
#include "../inc/fields.h"
#include "../inc/llg.h"
#include "../inc/sim.h"
#include "../inc/sf.h"
void sim::ramp_field(int argc,char *argv[])
{
    double rt=0.0,ef[3]={0,0,0};
    //rts = run timesteps
    unsigned int rts=0;
    std::string filestr;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Field details***" << std::endl;
    //The k-vector information is read in the file sf.cpp and sf_glob.cpp
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

    libconfig::Setting &setting = config::cfg.lookup("field");
    bool errstatus=setting.lookupValue("RunTime",rt);
    if(errstatus)
    {
        FIXOUT(config::Info,"Run time:" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run:" << rts << std::endl);
    }
    errstatus=setting.lookupValue("Temperature",llg::T);
    if(errstatus)
    {
        FIXOUT(config::Info,"Temperature:" << llg::T << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read temperature for sf (sf:Temperature (double))");
    }
    ef[0]=setting["Field"][0];
    ef[1]=setting["Field"][1];
    ef[2]=setting["Field"][2];
    //need to convert df to T per step
    errstatus=setting.lookupValue("OutputFilename",filestr);
    if(errstatus)
    {
        FIXOUT(config::Info,"Name of output file:" << filestr << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the name of the output file string (ramp_field:OutputFilename (string))");
    }
    std::ofstream ofs(filestr.c_str());
    if(!ofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file stream for writing");
    }
    llg::applied[0]=ef[0];
    llg::applied[1]=ef[1];
    llg::applied[2]=ef[2];
    for(unsigned int t = 0 ; t < rts ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
            ofs << static_cast<double>(t)*llg::dt << "\t" << llg::applied[0] << "\t" << llg::applied[1] << "\t" << llg::applied[2];
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                ofs << "\t" << spins::mag(i,0) << "\t" << spins::mag(i,1) << "\t" << spins::mag(i,2);
            }
            ofs << std::endl;
        }
        llg::integrate(t);
    }



}
