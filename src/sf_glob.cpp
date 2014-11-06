// File: sf_glob.cpp
// Note: originall dsf_glob.cpp
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 06 Nov 2014 12:01:47
#include "../inc/llg.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/sf.h"
#include "../inc/unitcell.h"
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

//Reads the parameters and sets up parts of the SF calculation
namespace sf
{
    //calculate the sf?
    bool csf=false;
    //The number of K-points we are
    //interested in
    unsigned int nk=0;
    //Ry and Rz are 3D rotation matrices and uo is the unitary operation
    Array2D<double> Ry,Rz,uo;
    Array2D<unsigned int> kpoints;
    Array<double> qa;
    //sfupdate=the sample timesteps in units of llg::update.
    unsigned int sfupdate=0;

    void readSFParam(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Dynamic structure factor details***" << std::endl;
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
        bool errstatus=false;
        std::string sffile;
        libconfig::Setting &setting = config::cfg.lookup("sf");
        errstatus=setting.lookupValue("Calculate",csf);
        //do some error handling
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read whether SF calculation is selected.");
        }
        else
        {
            FIXOUT(config::Info,"Outputting SF?:" << config::isTF(csf) << std::endl);
            if(csf==false)
            {
                //we do not need to any more processing of the dsf information
                return;
            }
            errstatus=setting.lookupValue("InputFile",sffile);
            if(errstatus==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You specified the use of the SF calculation, however you have to speficy an input file (sf:InputFile) in the config file.");
            }
            else
            {
                FIXOUT(config::Info,"SF input file:" << sffile << std::endl);
            }
        }
        try
        {
            config::cfg.readFile(sffile.c_str());
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
        libconfig::Setting &nsetting = config::cfg.lookup("sf");
        //if we have not exited or returned by this point we need to read the input file
        //First read the (predicted) quantization axis
        qa.resize(3);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            //initialize to zero
            qa[i]=0.0;
            //read from the config file
            qa[i]=nsetting["QuantizationAxis"][i];
        }
        FIXOUTVEC(config::Info,"Axis to calculate deviations around:",qa[0],qa[1],qa[2]);
        //resize the rotation matrices
        Ry.resize(3,3);Rz.resize(3,3);Ry.IFill(0);Rz.IFill(0);
        //calculate the theta (from z) and phi (from x) of the QA
        const double theta=acos(qa[2]);
        FIXOUT(config::Info,"Angle of quantization axis from z:" << theta*180.0/M_PI << " [degrees]" << std::endl);
        const double phi=atan2(qa[1],qa[0]);
        FIXOUT(config::Info,"Angle of quantization axis from x:" << phi*180.0/M_PI << " [degrees]" << std::endl);
        //Note here out rotation is towards the axis so
        //we take the negative
        Ry(0,0)=cos(-theta);Ry(0,1)=0;Ry(0,2)=sin(-theta);
        Ry(1,0)=0;Ry(1,1)=1;Ry(1,2)=0;
        Ry(2,0)=-sin(-theta);Ry(2,1)=0;Ry(2,2)=cos(-theta);
        Rz(0,0)=cos(-phi);Rz(0,1)=-sin(-phi);Rz(0,2)=0;
        Rz(1,0)=sin(-phi);Rz(1,1)=cos(-phi);Rz(1,2)=0;
        Rz(2,0)=0;Rz(2,1)=0;Rz(2,2)=1;
        config::printline(config::Info);
        FIXOUTVEC(config::Info,"Rotation matrix (Ry):",Ry(0,0),Ry(0,1),Ry(0,2));
        FIXOUTVEC(config::Info,"",Ry(1,0),Ry(1,1),Ry(1,2));
        FIXOUTVEC(config::Info,"",Ry(2,0),Ry(2,1),Ry(2,2));
        config::printline(config::Info);
        FIXOUTVEC(config::Info,"Rotation matrix (Rz):",Rz(0,0),Rz(0,1),Rz(0,2));
        FIXOUTVEC(config::Info,"",Rz(1,0),Rz(1,1),Rz(1,2));
        FIXOUTVEC(config::Info,"",Rz(2,0),Rz(2,1),Rz(2,2));
        config::printline(config::Info);

        errstatus=nsetting.lookupValue("NumberKPoints",nk);
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the number of K-points (sf:NumberKPoints (integer))");
        }
        FIXOUT(config::Info,"Number of K-points:" << nk << std::endl);
        kpoints.resize(nk,3);
        for(unsigned int i = 0 ; i < nk ; i++)
        {
            std::stringstream sstr;
            sstr << "KPoint" << i;
            std::string str=sstr.str();
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                kpoints(i,j)=nsetting[str.c_str()][j];
            }
            std::stringstream sstr1;sstr1 << "K-point " << i;
            str=sstr1.str();
            FIXOUTVEC(config::Info,str,kpoints(i,0),kpoints(i,1),kpoints(i,2));
        }
        uo.resize(geom::ucm.GetNMS(),3);
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            std::stringstream sstr;
            sstr << "UnitaryOperator" << i;
            std::string str=sstr.str();
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                uo(i,j)=nsetting[str.c_str()][j];
            }
            std::stringstream sstr1;
            sstr1 << "Unitary operator for magnetic species " << i << ":";
            str=sstr1.str();
            FIXOUTVEC(config::Info,str,uo(i,0),uo(i,1),uo(i,2));
            std::stringstream sstr2;
            if(uo(i,0)>0)
            {
                sstr2 << "Sx -> Sx";
            }
            else
            {
                sstr2 << "Sx -> -Sx";
            }
            if(uo(i,1)>0)
            {
                sstr2 << " , Sy -> Sy";
            }
            else
            {
                sstr2 << " , Sy -> -Sy";
            }
            if(uo(i,2)>0)
            {
                sstr2 << " , Sz -> Sz";
            }
            else
            {
                sstr2 << " , Sz -> -Sz";
            }
            str=sstr2.str();
            FIXOUT(config::Info,"" << str << std::endl);
        }
        errstatus=nsetting.lookupValue("Update",sfupdate);
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the sample time (sf:update) from the sf config file.");
        }
        else
        {
            FIXOUT(config::Info,"Sample timesteps (in units of spins:update):" << sfupdate << std::endl);
        }
    }
}
