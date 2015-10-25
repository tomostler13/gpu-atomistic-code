// File: sf_glob.cpp
// Note: originall dsf_glob.cpp
// Author:Tom Ostler
// Created: 23 Oct 2015
// Last-modified: 25 Oct 2015 21:44:57
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
    //update (in units of dipole field)
    unsigned int upd=0;
    //count for update
    unsigned int count=0;
    bool cab[3][3];
    Array2D<fftw_plan> rscfP;
    fftw_plan SxP,SyP,SzP;
    unsigned int nzpcplxdim=0;
    Array2D<unsigned int> rpoints;
    bool inc[3]={false,false,false};
    std::ofstream Coxx,Coxy,Coxz,Coyx,Coyy,Coyz,Cozx,Cozy,Cozz;
    //normalisation factor
    double normfac=1.0;
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
            nzpcplxdim=(geom::dim[2]*geom::Nk[2]/2)+1;
            normfac=(1./(static_cast<double>(geom::dim[0]*geom::dim[1]*geom::dim[2]*geom::Nk[0]*geom::Nk[1]*geom::Nk[2])))*(1./(static_cast<double>(geom::dim[0]*geom::dim[1]*geom::dim[2]*geom::Nk[0]*geom::Nk[1]*geom::Nk[2])));
            FIXOUT(config::Info,"Normalisation for correlation function:" << 1/normfac << std::endl);
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

            FIXOUT(config::Info,"Planning forward transforms of spin arrays for correlation function:" << std::flush);
            if(cab[0][0]==true || cab[0][1]==true || cab[0][2]==true || cab[1][0]==true || cab[2][0]==true)
            {
                inc[0]=true;
                spins::Srx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
                spins::Skx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                SxP=fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Srx.ptr(),spins::Skx.ptr(),FFTW_ESTIMATE);
            }
            if(cab[1][0]==true || cab[1][1]==true || cab[1][2]==true || cab[0][1]==true || cab[2][1]==true)
            {
                inc[1]=true;
                spins::Sky.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Sry.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
                SyP=fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Sry.ptr(),spins::Sky.ptr(),FFTW_ESTIMATE);
            }
            if(cab[2][0]==true || cab[2][1]==true || cab[2][2]==true || cab[0][2]==true || cab[1][2]==true)
            {
                inc[2]=true;
                spins::Skz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Srz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
                SzP=fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Srz.ptr(),spins::Skz.ptr(),FFTW_ESTIMATE);
            }
            config::Info << "Done" << std::endl;
            FIXOUT(config::Info,"Resizing Cr arrays and planning back transform of Sq . Sq*:" << std::flush);
            if(cab[0][0])
            {
                //do an in-place transform
                spins::Ckxx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Crxx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(0,0)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckxx.ptr(),spins::Crxx.ptr(),FFTW_ESTIMATE);
                Coxx.open("Cxx.dat");
                if(!Coxx.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Cxx correlation function.");
                }
            }
            if(cab[0][1])
            {
                spins::Ckxy.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Crxy.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(0,1)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckxy.ptr(),spins::Crxy.ptr(),FFTW_ESTIMATE);
                Coxy.open("Cxy.dat");
                if(!Coxy.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Cxy correlation function.");
                }
            }
            if(cab[0][2])
            {
                spins::Ckxz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Crxz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(0,2)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckxz.ptr(),spins::Crxz.ptr(),FFTW_ESTIMATE);
                Coxz.open("Cxz.dat");
                if(!Coxz.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Cxz correlation function.");
                }
            }
            if(cab[1][0])
            {
                spins::Ckyx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Cryx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(1,0)=rscfP(1,0)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckyx.ptr(),spins::Cryx.ptr(),FFTW_ESTIMATE);
                Coyx.open("Cyx.dat");
                if(!Coyx.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Cyx correlation function.");
                }
            }
            if(cab[1][1])
            {
                spins::Ckyy.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Cryy.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(1,1)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckyy.ptr(),spins::Cryy.ptr(),FFTW_ESTIMATE);
                Coyy.open("Cyy.dat");
                if(!Coyy.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Cyy correlation function.");
                }
            }
            if(cab[1][2])
            {
                spins::Ckyz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Cryz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(1,2)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckyz.ptr(),spins::Cryz.ptr(),FFTW_ESTIMATE);
                Coyz.open("Cyz.dat");
                if(!Coyz.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Cyz correlation function.");
                }
            }
            if(cab[2][0])
            {
                spins::Ckzx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Crzx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(2,0)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckzx.ptr(),spins::Crzx.ptr(),FFTW_ESTIMATE);
                Cozx.open("Czx.dat");
                if(!Cozx.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Czx correlation function.");
                }
            }
            if(cab[2][1])
            {
                spins::Ckzy.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Crzy.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(2,1)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckzy.ptr(),spins::Crzy.ptr(),FFTW_ESTIMATE);
                Cozy.open("Czy.dat");
                if(!Cozy.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Czy correlation function.");
                }
            }
            if(cab[2][2])
            {
                spins::Ckzz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
                spins::Crzz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

                rscfP(2,2)=fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],spins::Ckzz.ptr(),spins::Crzz.ptr(),FFTW_ESTIMATE);
                Cozz.open("Czz.dat");
                if(!Cozz.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open the file for writing the Czz correlation function.");
                }
            }
            config::Info << "Done" << std::endl;
            if(nsetting.lookupValue("NumberPoints",nr))
            {
                FIXOUT(config::Info,"Number of points to output correlation function for:" << nr << std::endl);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not read the number of point in real space");
            }
            if(nr > 0)
            {
                rpoints.resize(nr,3);
                for(unsigned int i = 0 ; i < nr ; i++)
                {
                    std::stringstream rpsstr;
                    rpsstr << "Ri_Minus_Rj_" << i;
                    std::string rpstr=rpsstr.str();
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {

                        try
                        {
                            rpoints(i,j)=nsetting[rpstr.c_str()][j];
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
                    FIXOUTVEC(config::Info,rpstr,rpoints(i,0),rpoints(i,1),rpoints(i,2));
                }
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("The number of points in real space for the correlation should be greater than 0");
            }
            if(nsetting.lookupValue("Update",upd))
            {
                FIXOUT(config::Info,"Update (in units of dipole field update):" << upd << std::endl);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not read the update of the correlation function");
            }
        }
    }

    void calcRSCF(unsigned int t)
    {
        count++;
        if(count%upd==0)
        {
            //first copy the spin info to the real space arrays
            if(inc[0])
            {
                spins::Skx.IFill(0);
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    unsigned int xyz[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                    spins::Srx(xyz[0],xyz[1],xyz[2])=spins::Sx[i];
                }
                fftw_execute(SxP);
            }
            if(inc[1])
            {
                spins::Sky.IFill(0);
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    unsigned int xyz[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                    spins::Sry(xyz[0],xyz[1],xyz[2])=spins::Sy[i];
                }
                fftw_execute(SyP);
            }
            if(inc[2])
            {
                spins::Skz.IFill(0);
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    unsigned int xyz[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                    spins::Srz(xyz[0],xyz[1],xyz[2])=spins::Sz[i];
                }
                fftw_execute(SzP);
            }
            //Perform the Sq_alpha . Sq_beta*
            if(cab[0][0])
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckxx(i,j,k)[0]=spins::Skx(i,j,k)[0]*spins::Skx(i,j,k)[0]+spins::Skx(i,j,k)[1]*spins::Skx(i,j,k)[1];
                            //Complex part
                            spins::Ckxx(i,j,k)[1]=spins::Skx(i,j,k)[1]*spins::Skx(i,j,k)[0]-spins::Skx(i,j,k)[0]*spins::Skx(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(0,0));
                //output the correlation function
                outputRSCF(Coxx,spins::Crxx,t);
            }
            if(cab[0][1])//Cxy
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckxy(i,j,k)[0]=spins::Skx(i,j,k)[0]*spins::Sky(i,j,k)[0]+spins::Skx(i,j,k)[1]*spins::Sky(i,j,k)[1];
                            //Complex part
                            spins::Ckxy(i,j,k)[1]=spins::Skx(i,j,k)[1]*spins::Sky(i,j,k)[0]-spins::Skx(i,j,k)[0]*spins::Sky(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(0,1));
                //output the correlation function
                outputRSCF(Coxy,spins::Crxy,t);
            }
            if(cab[0][2])//Cxz
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckxz(i,j,k)[0]=spins::Skx(i,j,k)[0]*spins::Skz(i,j,k)[0]+spins::Skx(i,j,k)[1]*spins::Skz(i,j,k)[1];
                            //Complex part
                            spins::Ckxz(i,j,k)[1]=spins::Skx(i,j,k)[1]*spins::Skz(i,j,k)[0]-spins::Skx(i,j,k)[0]*spins::Skz(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(0,2));
                //output the correlation function
                outputRSCF(Coxz,spins::Crxz,t);
            }
            if(cab[1][0])//Cyx
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckyx(i,j,k)[0]=spins::Sky(i,j,k)[0]*spins::Skx(i,j,k)[0]+spins::Sky(i,j,k)[1]*spins::Skx(i,j,k)[1];
                            //Complex part
                            spins::Ckyx(i,j,k)[1]=spins::Sky(i,j,k)[1]*spins::Skx(i,j,k)[0]-spins::Sky(i,j,k)[0]*spins::Skx(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(1,0));
                //output the correlation function
                outputRSCF(Coyx,spins::Cryx,t);
            }
            if(cab[1][1])//Cyy
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckyy(i,j,k)[0]=spins::Sky(i,j,k)[0]*spins::Sky(i,j,k)[0]+spins::Sky(i,j,k)[1]*spins::Sky(i,j,k)[1];
                            //Complex part
                            spins::Ckyy(i,j,k)[1]=spins::Sky(i,j,k)[1]*spins::Sky(i,j,k)[0]-spins::Sky(i,j,k)[0]*spins::Sky(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(1,1));
                //output the correlation function
                outputRSCF(Coyy,spins::Cryy,t);
            }
            if(cab[1][2])//Cyz
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckyz(i,j,k)[0]=spins::Sky(i,j,k)[0]*spins::Skz(i,j,k)[0]+spins::Sky(i,j,k)[1]*spins::Skz(i,j,k)[1];
                            //Complex part
                            spins::Ckyz(i,j,k)[1]=spins::Sky(i,j,k)[1]*spins::Skz(i,j,k)[0]-spins::Sky(i,j,k)[0]*spins::Skz(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(1,2));
                //output the correlation function
                outputRSCF(Coyz,spins::Cryz,t);
            }
            if(cab[2][0])//Czx
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckzx(i,j,k)[0]=spins::Skz(i,j,k)[0]*spins::Skx(i,j,k)[0]+spins::Skz(i,j,k)[1]*spins::Skx(i,j,k)[1];
                            //Complex part
                            spins::Ckzx(i,j,k)[1]=spins::Skz(i,j,k)[1]*spins::Skx(i,j,k)[0]-spins::Skz(i,j,k)[0]*spins::Skx(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(2,0));
                //output the correlation function
                outputRSCF(Cozx,spins::Crzx,t);
            }
            if(cab[2][1])//Czy
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckzy(i,j,k)[0]=spins::Skz(i,j,k)[0]*spins::Sky(i,j,k)[0]+spins::Skz(i,j,k)[1]*spins::Sky(i,j,k)[1];
                            //Complex part
                            spins::Ckzy(i,j,k)[1]=spins::Skz(i,j,k)[1]*spins::Sky(i,j,k)[0]-spins::Skz(i,j,k)[0]*spins::Sky(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(2,1));
                //output the correlation function
                outputRSCF(Cozy,spins::Crzy,t);
            }
            if(cab[2][2])//Czz
            {
                for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                        {
                            //S_{\alpha} \cdot S_{\beta}^{*} = (S_{\alpha}^r+i S_{\alpha}^c) \cdot (S_{\beta}^r-i S_{\beta}^c) = S_{\alpha}^r S_{\beta}^r + S_{\alpha}^c S_{\beta}^c + i(S_{\alpha}^c S_{\beta}^r - S_{\alpha}^r S_{\beta}^c)
                            //Real part
                            spins::Ckzz(i,j,k)[0]=spins::Skz(i,j,k)[0]*spins::Skz(i,j,k)[0]+spins::Skz(i,j,k)[1]*spins::Skz(i,j,k)[1];
                            //Complex part
                            spins::Ckzz(i,j,k)[1]=spins::Skz(i,j,k)[1]*spins::Skz(i,j,k)[0]-spins::Skz(i,j,k)[0]*spins::Skz(i,j,k)[1];
                        }
                    }
                }
                //perform back transform
                fftw_execute(rscfP(2,2));
                //output the correlation function
                outputRSCF(Cozz,spins::Crzz,t);
            }


            count=0;
        }
    }
    std::string conv(unsigned int c)
    {
        if(c==0)
        {
            return("x");
        }
        else if(c==1)
        {
            return("y");
        }
        else if(c==2)
        {
            return("z");
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Cartesian component not recognised. Should be equal to 0,1 or 2");
            //This should never be reached as error handling is dealt with within the function
            return(0);
        }
    }
    void outputRSCF(std::ofstream& ops,Array3D<double>& C,unsigned int t)
    {
        for(unsigned int i = 0 ; i < nr ; i++)
        {
            unsigned int xyz[3]={rpoints(i,0),rpoints(i,1),rpoints(i,2)};
            ops << t << "\t" << i << "\t" << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << "\t" << C(xyz[0],xyz[1],xyz[2])*normfac << std::endl;
        }
        ops << std::endl << std::endl;
    }
}
