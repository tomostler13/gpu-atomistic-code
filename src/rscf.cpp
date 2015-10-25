// File: sf_glob.cpp
// Note: originall dsf_glob.cpp
// Author:Tom Ostler
// Created: 23 Oct 2015
// Last-modified: 25 Oct 2015 18:06:19
#include "../inc/llg.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include "../inc/spins.h"
#include "../inc/rscf.h"
#include <sstream>
#include <fftw3.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

//Reads the parameters and sets up parts of the SF calculation
namespace rscf
{
    //calculate correlation function?
    bool ccf=false;
    //number of points
    unsigned int nr=0;
    bool cab[3][3];
    Array2D<fftw_plan> rscfP;
    fftw_plan SxP,SyP,SzP;
    Array2D<unsigned int> rpoints;
    void initRSCF(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Correlation function details***" << std::endl;
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
        std::string cffile;
        libconfig::Setting &setting = config::cfg.lookup("rscf");
        errstatus=setting.lookupValue("Calculate",ccf);
        //do some error handling
        if(errstatus==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read whether the correlation function calculation is selected.");
        }
        else
        {
            FIXOUT(config::Info,"Outputting correlation function?:" << config::isTF(ccf) << std::endl);
            if(ccf==false)
            {
                //we do not need to any more processing of the dsf information
                return;
            }
            errstatus=setting.lookupValue("InputFile",cffile);
            if(errstatus==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You specified the use of the real space correlation function calculation, however you have to speficy an input file (rscf:InputFile) in the config file.");
            }
            else
            {
                FIXOUT(config::Info,"Correlation function input file:" << cffile << std::endl);
            }
        }
        if(ccf)
        {
            try
            {
                config::cfg.readFile(cffile.c_str());
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
            libconfig::Setting &nsetting = config::cfg.lookup("rscf");
            rscfP.resize(3,3);
            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                try
                {
                    cab[0][i]=nsetting["Cx_beta"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Setting " << snf.getPath();
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
                try
                {
                    cab[1][i]=nsetting["Cy_beta"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Setting " << snf.getPath();
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
                try
                {
                    cab[2][i]=nsetting["Cz_beta"][i];
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
            FIXOUTVEC(config::Info,"Elements of rscf to be calculated:",config::isTF(cab[0][0]),config::isTF(cab[0][1]),config::isTF(cab[0][2]));
            FIXOUTVEC(config::Info,"",config::isTF(cab[1][0]),config::isTF(cab[1][1]),config::isTF(cab[1][2]));
            FIXOUTVEC(config::Info,"",config::isTF(cab[2][0]),config::isTF(cab[2][1]),config::isTF(cab[2][2]));


            //create the plans for the forward transform of the spins
            bool inc[3]={false,false,false};
            FIXOUT(config::Info,"Planning forward transforms of spin arrays for correlation function:" << std::flush);
            if(cab[0][0]==true || cab[0][1]==true || cab[0][2]==true || cab[1][0]==true || cab[2][0]==true)
            {
                inc[0]=true;
                spins::Srx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
                spins::Skx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                SxP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Srx.ptr(),spins::Skx.ptr(),FFTW_ESTIMATE);
            }
            if(cab[1][0]==true || cab[1][1]==true || cab[1][2]==true || cab[0][1]==true || cab[2][1]==true)
            {
                inc[1]=true;
                spins::Sky.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Sry.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
                SyP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Sry.ptr(),spins::Sky.ptr(),FFTW_ESTIMATE);
            }
            if(cab[2][0]==true || cab[2][1]==true || cab[2][2]==true || cab[0][2]==true || cab[1][2]==true)
            {
                inc[2]=true;
                spins::Skz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Srz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
                SzP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Srz.ptr(),spins::Skz.ptr(),FFTW_ESTIMATE);
            }
            config::Info << "Done" << std::endl;
            FIXOUT(config::Info,"Resizing Cr arrays and planning back transform of Sq . Sq*:" << std::flush);
            if(cab[0][0])
            {
                //do an in-place transform
                spins::Ckxx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Crxx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(0,0)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckxx.ptr(),spins::Crxx.ptr(),FFTW_ESTIMATE);
            }
            if(cab[0][1])
            {
                spins::Ckxy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Crxy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(0,1)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckxy.ptr(),spins::Crxy.ptr(),FFTW_ESTIMATE);
            }
            if(cab[0][2])
            {
                spins::Ckxz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Crxz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(0,2)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckxz.ptr(),spins::Crxz.ptr(),FFTW_ESTIMATE);
            }
            if(cab[1][0])
            {
                spins::Ckyx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Cryx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(1,0)=rscfP(1,0)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckyx.ptr(),spins::Cryx.ptr(),FFTW_ESTIMATE);
            }
            if(cab[1][1])
            {
                spins::Ckyy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Cryy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(1,1)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckyy.ptr(),spins::Cryy.ptr(),FFTW_ESTIMATE);
            }
            if(cab[1][2])
            {
                spins::Ckyz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Cryz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(1,2)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckyz.ptr(),spins::Cryz.ptr(),FFTW_ESTIMATE);
            }
            if(cab[2][0])
            {
                spins::Ckzx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Crzx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(2,0)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckzx.ptr(),spins::Crzx.ptr(),FFTW_ESTIMATE);
            }
            if(cab[2][1])
            {
                spins::Ckzy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Crzy.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(2,1)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckzy.ptr(),spins::Crzy.ptr(),FFTW_ESTIMATE);
            }
            if(cab[2][2])
            {
                spins::Ckzz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
                spins::Crzz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);

                rscfP(2,2)=fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],spins::Ckzz.ptr(),spins::Crzz.ptr(),FFTW_ESTIMATE);
            }
            config::Info << "Done" << std::endl;
        }
    }
}
