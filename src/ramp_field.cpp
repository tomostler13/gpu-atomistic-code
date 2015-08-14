// File: ramp_field.cpp
// Author: Tom Ostler
// Created: 13 May 2015
// Last-modified: 14 Aug 2015 12:51:52
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
    double rt=0.0,et=0.0,ef[3]={0,0,0},rf[3]={0,0,0};
    //rf is a unit vector to specify the direction in which the field
    //is changed
    //rts = run timesteps
    unsigned int rts=0,ets=0,numint=0;
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
    bool errstatus=false;setting.lookupValue("RunTime",rt);
    if(errstatus)
    {
        FIXOUT(config::Info,"Run time:" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run:" << rts << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the run time. Setting field.RunTime (double) in seconds");
    }
    if(setting.lookupValue("EquilTime",et))
    {
        FIXOUT(config::Info,"Equilibration time:" << et << std::endl);
        ets = static_cast<unsigned int>(et/llg::dt+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting field.EquilTime (double)");
    }
    if(setting.lookupValue("NumRateIntervals",numint))
    {
        FIXOUT(config::Info,"Number of rate intervals:" << numint << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the setting field.NumRateIntervals (int)");
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
    for(unsigned int i = 0 ; i < 3 ; i++)
    {
        try
        {
            ef[i]=setting["EqField"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting "<< snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
    }
    FIXOUTVEC(config::Info,"Equilibration field:",ef[0],ef[1],ef[2]);
    for(unsigned int i = 0 ; i < 3 ; i++)
    {
        try
        {
            rf[i]=setting["FieldChangeUnit"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting "<< snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
    }
    if(sqrt(rf[0]*rf[0]+rf[1]*rf[1]+rf[2]*rf[2])<0.99 && sqrt(rf[0]*rf[0]+rf[1]*rf[1]+rf[2]*rf[2]) > 1.01)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("FieldChangeUnit should be a unit vector, check your config file.");
    }
    else
    {
        FIXOUTVEC(config::Info,"Unit vector for field change:",rf[0],rf[1],ef[2]);
    }
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
    Array<double> fields,rates;
    fields.resize(numint+1);
    rates.resize(numint);
    fields.IFill(0);
    rates.IFill(0);
    for(unsigned int i = 0 ; i < numint+1 ; i++)
    {
        try
        {
            fields[i]=setting["fields"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
    }
    llg::applied[0]=ef[0];
    llg::applied[1]=ef[1];
    llg::applied[2]=ef[2];
    Array<double> dF;
    dF.resize(numint);
    dF.IFill(0);
    for(unsigned int i = 0 ; i < numint ; i++)
    {
        std::stringstream sstr,sstr1;
        sstr << "Interval " << i << " has field (magnitude) range ";
        sstr1 << fields[i] << " -> " << fields[i+1] << " [T]";
        std::string str=sstr.str(),str1=sstr1.str();
        FIXOUT(config::Info,str << str1 << std::endl);
        try
        {
            rates[i]=setting["Rates"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
        std::stringstream sstrn;
        sstr << "Interval " << i << " has rate:";
        std::string strn=sstrn.str();
        FIXOUT(config::Info,strn << rates[i] << std::endl);
        dF[i]=rates[i]*llg::dt;
        FIXOUT(config::Info,"Field increase (per update):" << dF[i] << " [T]" << std::endl);

    }
    std::ofstream Ttime("Field_time.dat");
    for(unsigned int t = 0 ; t < ets ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);

        }
        llg::integrate(t);
    }
    unsigned int time_counter=0;
    if(rates[0]>0)//then we are increasing field
    {
        for(unsigned int i = 0 ; i < numint ; i++)
        {
            for(double F = fields[i] ; F < fields[i+1] ; F+=dF[i])
            {
                if(time_counter%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(time_counter);
                    ofs << static_cast<double>(time_counter)*llg::dt << "\t" << llg::applied[0] << "\t" << llg::applied[1] << "\t" << llg::applied[2];
                    for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
                    {
                        ofs << "\t" << spins::mag(i,0) << "\t" << spins::mag(i,1) << "\t" << spins::mag(i,2);
                    }
                    ofs << std::endl;
                    llg::applied[0]=F*rf[0];
                    llg::applied[1]=F*rf[1];
                    llg::applied[2]=F*rf[2];
                }

                llg::integrate(time_counter);
                time_counter++;
            }
        }
    }

    ofs.close();
    if(ofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close output file.");
    }

}
