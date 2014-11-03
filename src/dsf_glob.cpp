// File: dsf_glob.cpp
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 03 Nov 2014 11:32:10
#include "../inc/llg.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/dsf.h"
#include "../inc/unitcell.h"
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

//Reads the parameters and sets up parts of the DSF calculation
namespace dsf
{
    //calculate the dsf?
    bool cdsf=false;
    //The number of symmetry directions we are
    //interested int
    unsigned int nsd=0;
    //Ry and Rz are 3D rotation matrices and uo is the unitary operation
    Array2D<double> Ry,Rz,uo;
    Array2D<unsigned int> symdir;
    Array<double> qa;
    //dsfupdate=the sample timesteps in units of llg::update.
    //ets = equilibration time
    //rts = run timesteps
    unsigned int dsfupdate=0,ets=0,rts=0;
    double et=0.0,rt=0.0,T=0.0;
    void readDSFParam(int argc,char *argv[])
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
        std::string dsffile;
        libconfig::Setting &setting = config::cfg.lookup("dsf");
        errstatus=setting.lookupValue("Calculate",cdsf);
        //do some error handling
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read whether DSF calculation is selected.");
        }
        else
        {
            FIXOUT(config::Info,"Outputting DSF?:" << config::isTF(cdsf) << std::endl);
            if(cdsf==false)
            {
                //we do not need to any more processing of the dsf information
                return;
            }
            errstatus=setting.lookupValue("InputFile",dsffile);
            if(errstatus==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You specified the use of the DSF calculation, however you have to speficy an input file (dsf:InputFile) in the config file.");
            }
            else
            {
                FIXOUT(config::Info,"DSF input file:" << dsffile << std::endl);
            }
        }
        try
        {
            config::cfg.readFile(dsffile.c_str());
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
        libconfig::Setting &nsetting = config::cfg.lookup("dsf");
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

        errstatus=nsetting.lookupValue("NumberOfSymmetryDirections",nsd);
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the number of symmetry directions (dsf:nsd)");
        }
        FIXOUT(config::Info,"Number of symmetry directions:" << nsd << std::endl);
        symdir.resize(nsd,3);
        for(unsigned int i = 0 ; i < nsd ; i++)
        {
            std::stringstream sstr;
            sstr << "SymmetryDirection" << i;
            std::string str=sstr.str();
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                symdir(i,j)=nsetting[str.c_str()][j];
            }
            std::stringstream sstr1;sstr1 << "Symmetry direction " << i;
            str=sstr1.str();
            FIXOUTVEC(config::Info,str,symdir(i,0),symdir(i,1),symdir(i,2));
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
        errstatus=nsetting.lookupValue("Update",dsfupdate);
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the sample time (dsf:update) from the dsf config file.");
        }
        else
        {
            FIXOUT(config::Info,"Sample timesteps (in units of spins:update):" << dsfupdate << std::endl);
        }
        errstatus=nsetting.lookupValue("EquilibrationTime",et);
        if(errstatus)
        {
            FIXOUT(config::Info,"Equilibration time:" << et << " [s]" << std::endl);
            ets=static_cast<unsigned int>((et/llg::dt)+0.5);
            FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
        }
        errstatus=nsetting.lookupValue("RunTime",rt);
        if(errstatus)
        {
            FIXOUT(config::Info,"Run time:" << rt << " [s]" << std::endl);
            rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
            FIXOUT(config::Info,"Number of timesteps for run:" << rts << std::endl);
        }
        errstatus=nsetting.lookupValue("Temperature",T);
        if(errstatus)
        {
            FIXOUT(config::Info,"Temperature:" << T << " [K]" << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read temperature for dsf (dsf:Temperature (double))");
        }
    }
}
