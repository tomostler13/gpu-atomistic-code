// File: anis.cpp
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 24 Sep 2014 11:06:14
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/anis.h"
#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>
//#define FIXOUT(a,b) a.width(75);a << std::left << b;
//#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;
//#define SUCCESS(a) a << "Done" << std::endl;
#include "../inc/defines.h"

namespace anis
{
    std::string units;
    unsigned int nfou=0;
    Array<double> FirstOrderUniaxK;
    Array2D<double> FirstOrderUniaxDir;
    void initAnis(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Anisotropy details***" << std::endl;
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
        libconfig::Setting &setting = config::cfg.lookup("anisotropy");
        try
        {
        }
        catch(...)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting type exception");
        }
        setting.lookupValue("units",units);
        if(units==""){std::cerr << "Units not read correctly\t" << units << std::endl;}
        FIXOUT(config::Info,"Anisotropy units:" << units << std::endl);
        setting.lookupValue("NumberFirstOrderUniaxAnis",nfou);
        FIXOUT(config::Info,"Number of first order uniaxial anisotropy axes" << nfou << std::endl);
        FirstOrderUniaxK.resize(nfou);
        FirstOrderUniaxDir.resize(nfou,3);
        for(unsigned int i = 0 ; i < nfou ; i++)
        {
            FirstOrderUniaxK[i]=0.0;
            std::stringstream sstrK,sstrD;
            std::string strK,strD;
            sstrK << "FirstOrderUniaxK_" << i+1;
            sstrD << "FirstOrderUniaxDir_" << i+1;
            strK=sstrK.str();
            strD=sstrD.str();
            setting.lookupValue(strK.c_str(),FirstOrderUniaxK[i]);
            FirstOrderUniaxDir(i,0)=setting[strD.c_str()][0];
            FirstOrderUniaxDir(i,1)=setting[strD.c_str()][1];
            FirstOrderUniaxDir(i,2)=setting[strD.c_str()][2];
            std::stringstream outsstr;
            std::string outstr;
            outsstr << "Uniaxial anisotropy constant " << i+1 << ":";
            outstr=outsstr.str();
            FIXOUT(config::Info,outstr.c_str() << FirstOrderUniaxK[i] << " " << units << std::endl);
            FIXOUTVEC(config::Info,"Direction:",FirstOrderUniaxDir(i,0),FirstOrderUniaxDir(i,1),FirstOrderUniaxDir(i,2));
        }
        SUCCESS(config::Info);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            if(units=="mRy")
            {
                FirstOrderUniaxK(i)*=2.179872172e-18;
                //now to milli
                FirstOrderUniaxK(i)*=1.0e-3;
            }
            else if(units=="joules" || units=="Joules" || units=="J")
            {
                //do nothing
            }
            else if(units=="eV")
            {
                FirstOrderUniaxK(i)*=1.602176565e-19;
            }
            else if(units=="meV")
            {
                FirstOrderUniaxK(i)*=1.602176565e-19;
                FirstOrderUniaxK(i)*=1e-3;
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Anisotropy units not recognised");
            }
        }
        for(unsigned int i = 0 ; i < nfou ; i++)
        {
//            FirstOrderUniaxK[i]*=(2.0/(mat::mu*mat::muB));
        }


    }
}
