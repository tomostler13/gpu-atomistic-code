// File: anis.cpp
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 18 Mar 2014 11:43:02
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
#include <string>
//#define FIXOUT(a,b) a.width(75);a << std::left << b;
//#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;
//#define SUCCESS(a) a << "Done" << std::endl;
#include "../inc/defines.h"

namespace anis
{
    Array3D<double> dT;
    Array2D<double> uniaxial_unit;

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
        dT.resize(mat::nspec,3,3);
        uniaxial_unit.resize(mat::nspec,3);
        uniaxial_unit.IFill(0);
        setting.lookupValue("units",units);
        if(units==""){std::cerr << "Units not read correctly\t" << units << std::endl;}
        FIXOUT(config::Info,"Anisotropy units:" << units << std::endl);
        FIXOUT(config::Info,"Reading anisotropy tensor:" << std::flush);
        for(unsigned int n = 0 ; n < mat::nspec ; n++)
        {
            std::stringstream sstrx,sstry,sstrz;
            sstrx << "d" << n;sstry << "d" << n;sstrz << "d" << n;
            sstrx << "xb";sstry << "yb";sstrz << "zb";
            std::string strx=sstrx.str(),stry=sstry.str(),strz=sstrz.str();

            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                dT(n,0,i)=setting[strx.c_str()][i];
                dT(n,1,i)=setting[stry.c_str()][i];
                dT(n,2,i)=setting[strz.c_str()][i];
            }
        }
        SUCCESS(config::Info);
        for(unsigned int n = 0 ; n < mat::nspec ; n++)
        {
            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                for(unsigned int j = 0 ; j < 3 ; j++)
                {
                    if(units=="mRy")
                    {
                        dT(n,i,j)*=2.179872172e-18;
                        //now to milli
                        dT(n,i,j)*=1.0e-3;
                    }
                    else if(units=="joules" || units=="Joules" || units=="J")
                    {
                        //do nothing
                    }
                    else if(units=="eV")
                    {
                        dT(n,i,j)*=1.602176565e-19;
                    }
                    else if(units=="meV")
                    {
                        dT(n,i,j)*=1.602176565e-19;
                        dT(n,i,j)*=1e-3;
                    }
                    else
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("Anisotropy units not recognised");
                    }
                }
            }
            FIXOUTVEC(config::Info,"dx_{beta}:",dT(n,0,0),dT(n,0,1),dT(n,0,2));
            FIXOUTVEC(config::Info,"dy_{beta}:",dT(n,1,0),dT(n,1,1),dT(n,1,2));
            FIXOUTVEC(config::Info,"dz_{beta}:",dT(n,2,0),dT(n,2,1),dT(n,2,2));
            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                for(unsigned int j = 0 ; j < 3 ; j++)
                {
                    dT(n,i,j)*=(2.0/(mat::muB*mat::mustore[n]));
                }
            }
            //for the use of the CSR neighbourlist we are onl considering uniaxial anisotropy (just diagonals)
            const double moddT=sqrt(dT(n,0,0)*dT(n,0,0)+dT(n,1,1)*dT(n,1,1)+dT(n,2,2)*dT(n,2,2));
            if(moddT>1e-32)
            {
                uniaxial_unit(n,0)=dT(n,0,0)/moddT;
                uniaxial_unit(n,1)=dT(n,1,1)/moddT;
                uniaxial_unit(n,2)=dT(n,2,2)/moddT;
            }
            else
            {
                uniaxial_unit(n,0)=0.0;
                uniaxial_unit(n,1)=0.0;
                uniaxial_unit(n,2)=0.0;
            }

        }

        if(config::useintmat==true)
        {
            if(mat::nspec<2)
            {

                FIXOUT(config::Info,"Normalising anisotropy and adding to interaction matrix:" << std::flush);
                intmat::Nrxx(0,0,0)+=(2.0*dT(0,0,0)/(mat::mustore[0]*mat::muB));
                intmat::Nrxy(0,0,0)+=(2.0*dT(0,0,1)/(mat::mustore[0]*mat::muB));
                intmat::Nrxz(0,0,0)+=(2.0*dT(0,0,2)/(mat::mustore[0]*mat::muB));
                intmat::Nryx(0,0,0)+=(2.0*dT(0,1,0)/(mat::mustore[0]*mat::muB));
                intmat::Nryy(0,0,0)+=(2.0*dT(0,1,1)/(mat::mustore[0]*mat::muB));
                intmat::Nryz(0,0,0)+=(2.0*dT(0,1,2)/(mat::mustore[0]*mat::muB));
                intmat::Nrzx(0,0,0)+=(2.0*dT(0,2,0)/(mat::mustore[0]*mat::muB));
                intmat::Nrzy(0,0,0)+=(2.0*dT(0,2,1)/(mat::mustore[0]*mat::muB));
                intmat::Nrzz(0,0,0)+=(2.0*dT(0,2,2)/(mat::mustore[0]*mat::muB));
                SUCCESS(config::Info);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You cannot currently use an interaction matrix with more than 1 species (Mar 2014)");
            }
        }


    }
}
