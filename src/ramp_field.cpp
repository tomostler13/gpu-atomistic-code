// File: ramp_field.cpp
// Author: Tom Ostler
// Created: 13 May 2015
// Last-modified: 10 Sep 2016 19:16:27
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
    double et=0.0,ef[3]={0,0,0},rf[3]={0,0,0};
    //rf is a unit vector to specify the direction in which the field
    //is changed
    unsigned int ets=0,numint=0;
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
    bool errstatus=false;
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
    bool ofsc=false;
    if(!setting.lookupValue("OutputFinalSpinConfig",ofsc))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read whether you wanted to output the final spin config. field.OutputFinalSpinConfig (bool)");
    }
    FIXOUT(config::Info,"Output final spin configuration?:" << config::isTF(ofsc) << std::endl);
    if(ofsc==true && geom::dim[0]*geom::dim[1]*geom::dim[2]>1)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The ouputting of the final spin config requires that the entire system is specified as a supercell (i.e. dim[3]={1,1,1};).");
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
            fields[i]=setting["Fields"][i];
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
    if(rates[0]<0)//then we are decreasing field
    {
        for(unsigned int i = 0 ; i < numint ; i++)
        {
            for(double F = fields[i] ; F > fields[i+1] ; F+=dF[i])
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
    util::calc_mag();
    util::output_mag(time_counter);
    if(ofsc)
    {
        std::ofstream finalucf("final_config.ucf");
        if(!finalucf.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not open file for outputting the final config file. The code will now attempt to output a vtu file that contains the information. This will have to be post processed.");
            unsigned int temp=-1;
            util::outputSpinsVTU(temp);
        }
        finalucf << geom::ucm.GetNMS() << std::endl;
        finalucf << geom::ucm.NumAtomsUnitCell() << std::endl;
        for(unsigned int i = 0 ; i < geom::ucm.NumAtomsUnitCell() ; i++)
        {
            finalucf << geom::ucm.GetSublattice(i) << "\t";
            for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
            {
                finalucf << geom::ucm.GetCoord(i,xyz)/static_cast<double>(geom::Nk[xyz]) << "\t";
            }
            finalucf << geom::ucm.GetMu(i) << "\t" << geom::ucm.GetDamping(i) << "\t" << geom::ucm.GetGamma(i) << "\t" << geom::ucm.GetElement(i) << "\t" << spins::Sx[i] << "\t" << spins::Sy[i] << "\t" << spins::Sz[i] << "\t" << geom::ucm.GetK1U(i) << "\t" << geom::ucm.GetK1UDir(i,0) << "\t" << geom::ucm.GetK1UDir(i,1) << "\t" << geom::ucm.GetK1UDir(i,2) << std::endl;
        }
        finalucf.close();
        if(finalucf.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close final_config.ucf file");
        }
    }
    ofs.close();
    if(ofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close output file.");
    }

}
