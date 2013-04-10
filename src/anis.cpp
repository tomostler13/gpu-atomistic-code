// File: anis.cpp
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 10 Apr 2013 13:46:54
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/mat.h"
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
    Array2D<double> dT;
    double uniaxial_unit[3]={0,0,0};

    std::string units;
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
        dT.resize(3,3);
        setting.lookupValue("units",units);
        if(units==""){std::cerr << "Units not read correctly\t" << units << std::endl;}
        FIXOUT(config::Info,"Anisotropy units:" << units << std::endl);
        FIXOUT(config::Info,"Reading anisotropy tensor:" << std::flush);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            dT(0,i)=setting["dxb"][i];
            dT(1,i)=setting["dyb"][i];
            dT(2,i)=setting["dzb"][i];
        }
        SUCCESS(config::Info);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                if(units=="mRy")
                {
                    dT(i,j)*=2.179872172e-18;
                    //now to milli
                    dT(i,j)*=1.0e-3;
                }
                else if(units=="joules" || units=="Joules" || units=="J")
                {
                    //do nothing
                }
                else if(units=="eV")
                {
                    dT(i,j)*=1.602176565e-19;
                }
                else if(units=="meV")
                {
                    dT(i,j)*=1.602176565e-19;
                    dT(i,j)*=1e-3;
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Anisotropy units not recognised");
                }
            }
        }
        FIXOUTVEC(config::Info,"dx_{beta}:",dT(0,0),dT(0,1),dT(0,2));
        FIXOUTVEC(config::Info,"dy_{beta}:",dT(1,0),dT(1,1),dT(1,2));
        FIXOUTVEC(config::Info,"dz_{beta}:",dT(2,0),dT(2,1),dT(2,2));
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                dT(i,j)*=(2.0/(mat::mu*mat::muB));
            }
        }
        //for the use of the CSR neighbourlist we are onl considering uniaxial anisotropy (just diagonals)
        const double moddT=sqrt(dT(0,0)*dT(0,0)+dT(1,1)*dT(1,1)+dT(2,2)*dT(2,2));
        uniaxial_unit[0]=dT(0,0)/moddT;
        uniaxial_unit[1]=dT(1,1)/moddT;
        uniaxial_unit[2]=dT(2,2)/moddT;
        if(config::useintmat)
        {

            FIXOUT(config::Info,"Normalising anisotropy and adding to interaction matrix:" << std::flush);
            /*        intmat::Nrxx(0,0,0)+=(2.0*dT[0][0]/(mat::mu*mat::muB));
                      intmat::Nrxy(0,0,0)+=(2.0*dT[0][1]/(mat::mu*mat::muB));
                      intmat::Nrxz(0,0,0)+=(2.0*dT[0][2]/(mat::mu*mat::muB));
                      intmat::Nryx(0,0,0)+=(2.0*dT[1][0]/(mat::mu*mat::muB));
                      intmat::Nryy(0,0,0)+=(2.0*dT[1][1]/(mat::mu*mat::muB));
                      intmat::Nryz(0,0,0)+=(2.0*dT[1][2]/(mat::mu*mat::muB));
                      intmat::Nrzx(0,0,0)+=(2.0*dT[2][0]/(mat::mu*mat::muB));
                      intmat::Nrzy(0,0,0)+=(2.0*dT[2][1]/(mat::mu*mat::muB));*/
            intmat::Nrzz(0,0,0)+=(2.0*dT(2,2)/(mat::mu*mat::muB));
            SUCCESS(config::Info);
        }


    }
}
