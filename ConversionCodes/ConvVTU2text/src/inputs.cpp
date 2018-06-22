// File: inputs.cpp
// Author:Tom Ostler
// Created: 22 Jun 2018
// Last-modified: 22 Jun 2018 09:21:45

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
    std::ofstream Info("convVTU2text.info");

    unsigned int nspins=0,lt=0,ut=0,dt=0;
    unsigned int dimin[3]={0,0,0};
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
        libconfig::Setting &setting = cfg.lookup("convVTU2text");
        if(setting.lookupValue("NumberSpin",nspins))
        {
            FIXOUT(Info,"Number of spins:" << nspins << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not get total number of spins");
        }
        if(setting.lookupValue("InputFilenameBase",fb))
        {
            FIXOUT(Info,"Base for filenames:" << fb << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the filename base (convVTU2text:InputFilenameBase (string))");
        }
        if(setting.lookupValue("StartTime",lt))
        {
            FIXOUT(Info,"Start time index:" << lt << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read start time index (convVTU2text:StartTime (int))");
        }
        if(setting.lookupValue("EndTime",ut))
        {
            FIXOUT(Info,"End time index:" << ut << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read end time index (convVTU2text:EndTime (int))");
        }
        if(setting.lookupValue("dTime",dt))
        {
            FIXOUT(Info,"Step in time index:" << lt << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read step in time index (convVTU2text:dTime (int))");
        }


    }
}
