// File: ramp_field.cpp
// Author: Tom Ostler
// Created: 13 May 2015
// Last-modified: 13 May 2015 20:33:28
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
    double et=0.0,rt=0.0,ef[3]={0.0,0.0,0.0},df[3]={0.0,0.0,0.0};
    //ets = equilibration time
    //rts = run timesteps
    //ef = equilibration field
    //df = change in field (T/ns)
    unsigned int ets=0,rts=0;
    std::string filestr;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Field Ramp details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("ramp_field");
    bool errstatus=setting.lookupValue("EquilibrationTime",et);
    if(errstatus)
    {
        FIXOUT(config::Info,"Equilibration time:" << et << " [s]" << std::endl);
        ets=static_cast<unsigned int>((et/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    }
    errstatus=setting.lookupValue("RunTime",rt);
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
    ef[0]=setting["EquilibrationField"][0];
    ef[1]=setting["EquilibrationField"][1];
    ef[2]=setting["EquilibrationField"][2];
    FIXOUTVEC(config::Info,"Field during equilibration [T]:",ef[0],ef[1],ef[2]);
    df[0]=setting["FieldChange"][0];
    df[1]=setting["FieldChange"][1];
    df[2]=setting["FieldChange"][2];
    FIXOUTVEC(config::Info,"Rate of field change [T/ns]:",df[0],df[1],df[2]);
    //need to convert df to T per step
    df[0]=df[0]/(1e-9/llg::dt);
    df[1]=df[1]/(1e-9/llg::dt);
    df[2]=df[2]/(1e-9/llg::dt);
    FIXOUTVEC(config::Info,"New change in field [T/step]:",df[0],df[1],df[2]);
    filestr=setting.lookupValue("OutputFilename",filestr);
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
    for(unsigned int t = 0 ; t < ets ; t++)
    {
        if(t%spins::update==0)
        {
            spins::mag.IFill(0);
            //add up the magnetization for each sublattice
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                unsigned int sl = geom::sublattice[i];
                spins::mag(sl,0)+=spins::Sx[i];
                spins::mag(sl,1)+=spins::Sy[i];
                spins::mag(sl,2)+=spins::Sz[i];
            }
            //divide by the number of spins in each sublattice
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                double oones=1./static_cast<double>(geom::ucm.GetNES(i));
                spins::mag(i,0)*=oones;
                spins::mag(i,1)*=oones;
                spins::mag(i,2)*=oones;
            }

            ofs << static_cast<double>(t)*llg::dt << "\t" << llg::applied[0] << "\t" << llg::applied[1] << "\t" << llg::applied[2];
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                ofs << "\t" << spins::mag(i,0) << "\t" << spins::mag(i,1) << "\t" << spins::mag(i,2);
            }
            ofs << std::endl;
        }
        llg::integrate(t);
    }


    for(unsigned int t = ets ; t < (ets+rts) ; t++)
    {
        llg::applied[0]+=df[0];
        llg::applied[1]+=df[1];
        llg::applied[2]+=df[2];
        if(t%spins::update==0)
        {
            spins::mag.IFill(0);
            //add up the magnetization for each sublattice
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                unsigned int sl = geom::sublattice[i];
                spins::mag(sl,0)+=spins::Sx[i];
                spins::mag(sl,1)+=spins::Sy[i];
                spins::mag(sl,2)+=spins::Sz[i];
            }
            //divide by the number of spins in each sublattice
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                double oones=1./static_cast<double>(geom::ucm.GetNES(i));
                spins::mag(i,0)*=oones;
                spins::mag(i,1)*=oones;
                spins::mag(i,2)*=oones;
            }

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
