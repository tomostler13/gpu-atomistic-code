// File: anis.cpp
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 28 Apr 2023 04:15:39 PM
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
    Array<double> k1u;
    Array2D<double> k1udir;
    double k2par,k2perp,k4par,k4perp;
    Array<double> k2perpdir,k2pardir;
    Array2D<double> k4perpdirs,k4pardirs;

    void readGlobalAnis()
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Global Anisotropy details***" << std::endl;
        bool anisset=false;
        if(!config::cfg.exists("anis"))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Setting class anis does not exist, setting defaults (all off).");
            k2perp = 0.0;
            k2par = 0.0;
            k4par = 0.0;
            k4perp = 0.0;
            k2perpdir.IFill(0);
            k2pardir.IFill(0);
            k4perpdirs.IFill(0);
            k4pardirs.IFill(0);
            FIXOUT(config::Info,"K2_Perp:" << k2perp << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K2 perp direction:",k2perpdir[0],k2perpdir[1],k2perpdir[2]);
            FIXOUT(config::Info,"K2_Parallel:" << k2par << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K2 parallel direction:",k2pardir[0],k2pardir[1],k2pardir[2]);
            FIXOUT(config::Info,"K4_Perp:" << k4perp << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K4 perp direction 1:",k4perpdirs(0,0),k4perpdirs(0,1),k4perpdirs(0,2));
            FIXOUTVEC(config::Info,"K4 perp direction 1:",k4perpdirs(0,0),k4perpdirs(0,1),k4perpdirs(0,2));
            FIXOUTVEC(config::Info,"K4 perp direction 2:",k4perpdirs(1,0),k4perpdirs(1,1),k4perpdirs(1,2));
            FIXOUTVEC(config::Info,"K4 perp direction 3:",k4perpdirs(2,0),k4perpdirs(2,1),k4perpdirs(2,2));
            FIXOUT(config::Info,"K4_Parallel:" << k4perp << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K4 parallel direction 1:",k4pardirs(0,0),k4pardirs(0,1),k4pardirs(0,2));
            FIXOUTVEC(config::Info,"K4 parallel direction 2:",k4pardirs(1,0),k4pardirs(1,1),k4pardirs(1,2));
            FIXOUTVEC(config::Info,"K4 parallel direction 3:",k4pardirs(2,0),k4pardirs(2,1),k4pardirs(2,2));
        }
        else
        {
            //Need to run a check to make sure that all of the magnetic moments are the same
            //because this part of the code has been hard coded!!!!!
            const double basemm=geom::ucm.GetMu(0);
            bool allEqual=true;
            for(unsigned int t = 1 ; t < geom::ucm.NumAtomsUnitCell() ; t++)
            {
                if(fabs((geom::ucm.GetMu(t)-basemm)/basemm)>1e-6)
                {
                    allEqual=false;
                    break;
                }
            }
            if(!allEqual)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Because the K2 and K4 terms have been hard coded in a less general way, the magnetic moments of all magnetic species MUST be the same. This needs fixing in the future to make the code more general");
            }
            libconfig::Setting &setting = config::cfg.lookup("anis");
            if(!setting.lookupValue("k2perp",k2perp))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Setting anis set, but no value for k2perp set (double)");
            }
            if(!setting.lookupValue("k2par",k2par))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Setting anis set, but no value for k2par set (double)");
            }
            if(!setting.lookupValue("k4par",k4par))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Setting anis set, but no value for k4par set (double)");
            }
            if(!setting.lookupValue("k4perp",k4perp))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Setting anis set, but no value for k4perp set (double)");
            }
            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                try
                {
                    k2pardir[i] = setting["k2pardir"][i];
                }
                catch(const libconfig::SettingNotFoundException & snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Exchange striction not set up properly. Setting " << snf.getPath() << " must be set.";
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
                try
                {
                    k2perpdir[i] = setting["k2perpdir"][i];
                }
                catch(const libconfig::SettingNotFoundException & snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Exchange striction not set up properly. Setting " << snf.getPath() << " must be set.";
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }

                std::stringstream k4pardir_sstr;
                std::stringstream k4perpdir_sstr;
                k4pardir_sstr << "k4pardir" << i+1;
                k4perpdir_sstr << "k4perpdir" << i+1;
                std::string k4pardir_str=k4pardir_sstr.str();
                std::string k4perpdir_str=k4perpdir_sstr.str();
                for(unsigned int j = 0 ; j < 3 ; j++)
                {
                    try
                    {
                        k4pardirs(i,j) = setting[k4pardir_str][j];
                    }
                    catch(const libconfig::SettingNotFoundException & snf)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream errsstr;
                        errsstr << "Setting not found exception caught. Exchange striction not set up properly. Setting " << snf.getPath() << " must be set.";
                        std::string errstr=errsstr.str();
                        error::errMessage(errstr);
                    }
                    try
                    {
                        k4perpdirs(i,j) = setting[k4perpdir_str][j];
                    }
                    catch(const libconfig::SettingNotFoundException & snf)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream errsstr;
                        errsstr << "Setting not found exception caught. Exchange striction not set up properly. Setting " << snf.getPath() << " must be set.";
                        std::string errstr=errsstr.str();
                        error::errMessage(errstr);
                    }


                }
            }
            FIXOUT(config::Info,"K2_Perp:" << k2perp << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K2 perp direction:",k2perpdir[0],k2perpdir[1],k2perpdir[2]);
            FIXOUT(config::Info,"K2_Parallel:" << k2par << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K2 parallel direction:",k2pardir[0],k2pardir[1],k2pardir[2]);
            FIXOUT(config::Info,"K4_Perp:" << k4perp << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K4 perp direction 1:",k4perpdirs(0,0),k4perpdirs(0,1),k4perpdirs(0,2));
            FIXOUTVEC(config::Info,"K4 perp direction 2:",k4perpdirs(1,0),k4perpdirs(1,1),k4perpdirs(1,2));
            FIXOUTVEC(config::Info,"K4 perp direction 3:",k4perpdirs(2,0),k4perpdirs(2,1),k4perpdirs(2,2));
            FIXOUT(config::Info,"K4_Parallel:" << k4perp << " [J]" << std::endl);
            FIXOUTVEC(config::Info,"K4 parallel direction 1:",k4pardirs(0,0),k4pardirs(0,1),k4pardirs(0,2));
            FIXOUTVEC(config::Info,"K4 parallel direction 2:",k4pardirs(1,0),k4pardirs(1,1),k4pardirs(1,2));
            FIXOUTVEC(config::Info,"K4 parallel direction 3:",k4pardirs(2,0),k4pardirs(2,1),k4pardirs(2,2));

        }//end of else
    }//end of ReadGlobalAnis function
}
