// File: inputs.cpp
// Author:Tom Ostler
// Created: 25 Feb 2016
// Last-modified: 24 Mar 2017 15:16:51

//Reads the config file
#include <cmath>
#include <iostream>
#include <libconfig.h++>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "../../../inc/error.h"
#include "../../../inc/array.h"
#include "../../../inc/array2d.h"
#include "../../../inc/array3d.h"
#include "../../../inc/array4d.h"
#include "../../../inc/defines.h"
#include "../../../inc/random.h"
#include "../inc/inputs.h"
namespace inputs
{
    //for outputting info
    std::ofstream Info("rscf.info");

    unsigned int nspins=0,lt=0,ut=0,dt=0;
    unsigned int dimin[3]={0,0,0},dimout[3]={0,0,0},Nk[3]={0,0,0};
    double abc[3] = {0.0,0.0,0.0};
    std::string fb;
    void readcff(int argc,char *argv[])
    {
        if(!Info.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for outputting information.");
        }
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
        libconfig::Setting &setting = cfg.lookup("rscf");
        if(setting.lookupValue("NumberSpin",nspins))
        {
            FIXOUT(Info,"Number of spins (in):" << nspins << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not get total number of spins");
        }
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            try
            {
                dimin[i] = setting["NCellsIn"][i];
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
        FIXOUTVEC(Info,"Dimensions read in:",dimin[0],dimin[1],dimin[2]);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            try
            {
                abc[i] = setting["abc"][i];
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
        FIXOUTVEC(Info,"Size of unit cell (abc):",abc[0],abc[1],abc[2]);
        if(setting.lookupValue("InputFilenameBase",fb))
        {
            FIXOUT(Info,"Base for filenames:" << fb << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the filename base (rscf:InputFilenameBase (string))");
        }
        if(setting.lookupValue("StartTime",lt))
        {
            FIXOUT(Info,"Start time index:" << lt << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read start time index (rscf:StartTime (int))");
        }
        if(setting.lookupValue("EndTime",ut))
        {
            FIXOUT(Info,"End time index:" << ut << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read end time index (rscf:EndTime (int))");
        }
        if(setting.lookupValue("dTime",dt))
        {
            FIXOUT(Info,"Step in time index:" << lt << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read step in time index (rscf:dTime (int))");
        }


    }
}
