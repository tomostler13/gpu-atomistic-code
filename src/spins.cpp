// File: spins.cpp
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 08 Jan 2018 20:28:01
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/intmat.h"
#include "../inc/defines.h"
#include "../inc/maths.h"
#include "../inc/random.h"
namespace spins
{
    Array5D<fftw_complex> Sk;
    Array5D<fftw_complex> Sr;
    Array4D<fftw_complex> hSk,hSr;
    Array4D<fftw_complex> dipSr;
    Array4D<fftw_complex> dipSk;
    //Possible correlation functions
    Array3D<double> Srx,Sry,Srz;
    //Similarly for complex
    Array3D<fftw_complex> Skx,Sky,Skz;
    Array3D<fftw_complex> Ckxx,Ckxy,Ckxz,Ckyx,Ckyy,Ckyz,Ckzx,Ckzy,Ckzz,Cksds;
    Array3D<double> Crxx,Crxy,Crxz,Cryx,Cryy,Cryz,Crzx,Crzy,Crzz,Crsds;
    Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
    Array2D<double> mag;
    Array2D<double> magEachAtUC;
    fftw_plan SP,dSP;
    unsigned int update=0,mag_calc_method=0;
    bool output_mag=true,mapout=false;
    std::ifstream sfs;
    void initSpins()
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Spin details***" << std::endl;
        FIXOUT(config::Info,"Resizing arrays and initializing spin directions:" << std::flush);
        mag.resize(geom::ucm.GetNMS(),3);
        mag.IFill(0);
        Sx.resize(geom::nspins);
        Sy.resize(geom::nspins);
        Sz.resize(geom::nspins);
        Sx.IFill(geom::nspins);
        Sy.IFill(geom::nspins);
        Sz.IFill(geom::nspins);

        SUCCESS(config::Info);
        FIXOUT(config::Info,"Method for calculating the magnetization:" << mag_calc_method << " (see notes in src/util.cpp, function -> calc_mag)" << std::endl);
        FIXOUT(config::Info,"Output magnetization?" << config::isTF(output_mag) << std::endl);
        if(config::exchm==0)
        {
            Sk.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            Sr.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            Sr.IFill(0);
            Sk.IFill(0);

            FIXOUT(config::Info,"Planning forward and transform of spin map:" << std::flush);
            int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zps;
            int idist=geom::zps;

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of spins (forward transform)" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << geom::ucm.GetNMS()*3 << std::endl);
            FIXOUT(config::Log,"Pointer of real space spins (Sr):" << Sr.ptr() << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUT(config::Log,"Pointer of reciprocal space spins (Sk):" << Sk.ptr() << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) " << FFTW_FORWARD << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            SP = fftw_plan_many_dft(3,n,geom::ucm.GetNMS()*3,Sr.ptr(),inembed,istride,idist,Sk.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);
            //forward transform of spin arrays
            SUCCESS(config::Info);
        }
        else if(config::dipm==0 && config::inc_dip==true)
        {
            dipSk.resize(3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            dipSr.resize(3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            dipSr.IFill(0);
            dipSk.IFill(0);

            FIXOUT(config::Info,"Method for calculating the magnetization:" << mag_calc_method << " (see notes in src/util.cpp, function -> calc_mag)" << std::endl);
            FIXOUT(config::Info,"Planning forward and transform of spin map:" << std::flush);
            int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zps;
            int idist=geom::zps;

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of spins (forward transform)" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << 3 << std::endl);
            FIXOUT(config::Log,"Pointer of real space spins (dipSr):" << dipSr.ptr() << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUT(config::Log,"Pointer of reciprocal space spins (dipSk):" << dipSk.ptr() << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) " << FFTW_FORWARD << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            dSP = fftw_plan_many_dft(3,n,3,dipSr.ptr(),inembed,istride,idist,dipSk.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_PATIENT);
            //forward transform of spin arrays
            SUCCESS(config::Info);
        }

        if(config::exchm>98)
        {
            hSk.resize(3,geom::dim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            hSr.resize(3,geom::dim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            hSr.IFill(0);
            hSk.IFill(0);

            FIXOUT(config::Info,"Planning forward and transform of set of 2D spin maps (layer by layer):" << std::flush);
            int n[2]={geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zpdim[1]*geom::Nk[1]*geom::zpdim[2]*geom::Nk[2];
            int idist=odist;

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of spins (forward transform)" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],0);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 2 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << 3*geom::dim[0]*geom::Nk[0] << std::endl);
            FIXOUT(config::Log,"Pointer of real space spins (hSr):" << hSr.ptr() << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],0);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUT(config::Log,"Pointer of reciprocal space spins (hSk):" << hSk.ptr() << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],0);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) " << FFTW_FORWARD << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            SP = fftw_plan_many_dft(2,n,3*geom::dim[0]*geom::Nk[0],hSr.ptr(),inembed,istride,idist,hSk.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);
            //forward transform of spin arrays
            SUCCESS(config::Info);

        }

    }
    void FFTForward()
    {
        Sk.IFill(0);
        Sr.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            unsigned int sl=geom::sublattice[i];
            Sr(sl,0,lc[0],lc[1],lc[2])[0]=Sx[i];
            Sr(sl,1,lc[0],lc[1],lc[2])[0]=Sy[i];
            Sr(sl,2,lc[0],lc[1],lc[2])[0]=Sz[i];
        }
        fftw_execute(SP);
    }
    void hFFTForward()
    {
        hSk.IFill(0);
        hSr.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            hSr(0,lc[0],lc[1],lc[2])[0]=Sx[i];
            hSr(1,lc[0],lc[1],lc[2])[0]=Sy[i];
            hSr(2,lc[0],lc[1],lc[2])[0]=Sz[i];
        }
        fftw_execute(SP);
    }
    void dipFFTForward()
    {
        dipSk.IFill(0);
        dipSr.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            double relmu=geom::mu[i];
            dipSr(0,lc[0],lc[1],lc[2])[0]=Sx[i]*relmu;
            dipSr(1,lc[0],lc[1],lc[2])[0]=Sy[i]*relmu;
            dipSr(2,lc[0],lc[1],lc[2])[0]=Sz[i]*relmu;
        }
        fftw_execute(dSP);
    }
    void eFFTForward()
    {
        Sk.IFill(0);
        Sr.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            unsigned int sl=geom::sublattice[i];
            Sr(sl,0,lc[0],lc[1],lc[2])[0]=eSx[i];
            Sr(sl,1,lc[0],lc[1],lc[2])[0]=eSy[i];
            Sr(sl,2,lc[0],lc[1],lc[2])[0]=eSz[i];
        }
        fftw_execute(SP);
    }
    void heFFTForward()
    {
        hSk.IFill(0);
        hSr.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            hSr(0,lc[0],lc[1],lc[2])[0]=eSx[i];
            hSr(1,lc[0],lc[1],lc[2])[0]=eSy[i];
            hSr(2,lc[0],lc[1],lc[2])[0]=eSz[i];
        }
        fftw_execute(SP);
    }
    void dipeFFTForward()
    {
        dipSk.IFill(0);
        dipSr.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            double relmu=geom::mu[i];
            dipSr(0,lc[0],lc[1],lc[2])[0]=eSx[i]*relmu;
            dipSr(1,lc[0],lc[1],lc[2])[0]=eSy[i]*relmu;
            dipSr(2,lc[0],lc[1],lc[2])[0]=eSz[i]*relmu;
        }
        fftw_execute(dSP);
    }
    void setSpinsConfig()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            //get which atom in the unit cell spin i is
            unsigned int aiuc=geom::lu(i,4);
            Sx[i]=geom::ucm.GetInitS(aiuc,0);
            Sy[i]=geom::ucm.GetInitS(aiuc,1);
            Sz[i]=geom::ucm.GetInitS(aiuc,2);
        }
    }
    void setSpinsRandom()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            double v1=0,v2=0,s=2.0,ss=0.0;
            while(s>1.0)
            {
                v1=2.0*Random::rand()-1.0;
                v2=2.0*Random::rand()-1.0;
                s=v1*v1+v2*v2;
            }
            ss=sqrt(1.0-s);
            Sx[i]=2.0*v1*ss;
            Sy[i]=2.0*v2*ss;
            Sz[i]=1.0-2.0*s;
        }
    }
    void setSpinsChequerX()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            int mycoords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            if((mycoords[0]+mycoords[1]+mycoords[2])%2==0)
            {//then we are on sublattice A
                spins::Sx[i]=1.0;
                spins::Sy[i]=0.0;
                spins::Sz[i]=0.0;
            }
            else
            {//then we are on sublattice B
                spins::Sx[i]=-1.0;
                spins::Sy[i]=0.0;
                spins::Sz[i]=0.0;
            }
        }
    }
    void setSpinsChequerY()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            int mycoords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            if((mycoords[0]+mycoords[1]+mycoords[2])%2==0)
            {//then we are on sublattice A
                spins::Sx[i]=0.0;
                spins::Sy[i]=1.0;
                spins::Sz[i]=0.0;
            }
            else
            {//then we are on sublattice B
                spins::Sx[i]=0.0;
                spins::Sy[i]=-1.0;
                spins::Sz[i]=0.0;
            }
        }
    }
    void setSpinsChequerZ()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            int mycoords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            if((mycoords[0]+mycoords[1]+mycoords[2])%2==0)
            {//then we are on sublattice A
                spins::Sx[i]=0.0;
                spins::Sy[i]=0.0;
                spins::Sz[i]=1.0;
            }
            else
            {//then we are on sublattice B
                spins::Sx[i]=0.0;
                spins::Sy[i]=0.0;
                spins::Sz[i]=-1.0;
            }
        }
    }
    void setSpinsResText(libconfig::Setting& setting)
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            //initialise all spins to zero and at the end check that
            //their modulus is NOT 0
            spins::Sx[i]=0.0;
            spins::Sy[i]=0.0;
            spins::Sz[i]=0.0;
        }
        std::string fname;
        if(setting.lookupValue("SpinFile",fname))
        {
            FIXOUT(config::Info,"Initialising spins from file:" << fname << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the file (SpinFile) to initialise spins");
        }
        std::ifstream ifs(fname.c_str());
        if(ifs.is_open())
        {
            unsigned int tmpns=0;
            ifs >> tmpns;
            if(tmpns!=geom::nspins)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("The number of spins read in is not the same as the number you are now trying to create.");
            }
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                ifs >> spins::Sx[i] >> spins::Sy[i] >> spins::Sz[i];
                const double mods=sqrt(spins::Sx[i]*spins::Sx[i] + spins::Sy[i]*spins::Sy[i] + spins::Sz[i]*spins::Sz[i]);
                if(mods < 1e-12)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Read a spin with a zero (<1e-12) spin lenth)");
                }
                spins::Sx[i]/=mods;
                spins::Sy[i]/=mods;
                spins::Sz[i]/=mods;
            }
            ifs.close();
            if(ifs.is_open())
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not close spin file used for initialising spins.");
            }
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open spin file to initialise spin, check your filename.");
        }
    }

    void setSpinsSpecies(libconfig::Setting& setting)
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            //initialise all spins to zero and at the end check that
            //their modulus is NOT 0
            spins::Sx[i]=0.0;
            spins::Sy[i]=0.0;
            spins::Sz[i]=0.0;
        }
        std::string *initm=NULL;
        initm = new std::string [geom::ucm.GetNMS()];
        //std::cout << "Number of species = " << geom::ucm.GetNMS() << std::endl;
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            double chequergridscale[3]={0,0,0};
            std::stringstream sstr;
            sstr << "SublatticeInit" << i+1;
            std::string str=sstr.str();
            if(setting.lookupValue(str.c_str(),initm[i]))
            {
                if(initm[i]!="chequer" && initm[i]!="uniformz")
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "You have chosen a value for spin initialisation (read as " << initm[i] << " ) on a sublattice by sublattice basis that is not recognised. So far the only ones implemented are: chequer and uniformz";
                    std::string errstr=errsstr.str();

                    error::errMessage(errstr.c_str());
                }
                FIXOUT(config::Info,str << initm[i] << std::endl);
                if(initm[i]=="chequer")
                {
                    std::stringstream gsstr;
                    gsstr << "ChequerGridScale" << i+1;
                    std::string gstr=gsstr.str();
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        try
                        {
                            chequergridscale[xyz]=setting[gstr.c_str()][xyz];
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
                }
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                std::stringstream errsstr;
                errsstr << "Could not read the initialisation method for sublattice " << i+1 << ". llg.SublatticeInit (std::string)";
                std::string errstr=errsstr.str();
                error::errMessage(errstr.c_str());
            }
            for(unsigned int spin = 0 ; spin < geom::nspins ; spin++)
            {
                if(initm[i]=="uniformz" && geom::lu(spin,3)==static_cast<int>(i))
                {
                    spins::Sx[spin]=0;
                    spins::Sy[spin]=0;
                    spins::Sz[spin]=1.0;
//                    std::cout << "Spin " << spin << " is species " << i << std::endl;
                }
                else if(initm[i]=="chequer" && geom::lu(spin,3)==static_cast<int>(i))
                {
//                    std::cout << "Spin " << spin << " is species " << i << std::endl;
                    int mycoords[3]={static_cast<int>(static_cast<double>(geom::lu(spin,0))*chequergridscale[0]+0.5),static_cast<int>(static_cast<double>(geom::lu(spin,1))*chequergridscale[1]+0.5),static_cast<int>(static_cast<double>(geom::lu(spin,2))*chequergridscale[2]+0.5)};
                    if((mycoords[0]+mycoords[1]+mycoords[2])%2==0)
                    {//then we are on sublattice A
                    //std::cout << "Sublattice A " << mycoords[0] << "\t" << mycoords[1] << "\t" << mycoords[2] << std::endl;
                        spins::Sx[spin]=0.0;
                        spins::Sy[spin]=0.0;
                        spins::Sz[spin]=1.0;
                    }
                    else
                    {//then we are on sublattice B
                    //std::cout << "Sublattice B " << mycoords[0] << "\t" << mycoords[1] << "\t" << mycoords[2] << std::endl;
                        spins::Sx[spin]=0.0;
                        spins::Sy[spin]=0.0;
                        spins::Sz[spin]=-1.0;
                    }
                }
            }
        }
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            if(sqrt(spins::Sx[i]*spins::Sx[i]+spins::Sy[i]*spins::Sy[i]+spins::Sz[i]*spins::Sz[i])<1e-10)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Error in spin initialisation, found a spin with zero modulus (all should have a value of 1)");
            }
        }

    }
    void setSpinsVampire(libconfig::Setting& setting)
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            //initialise all spins to zero and at the end check that
            //their modulus is NOT 0
            spins::Sx[i]=0.0;
            spins::Sy[i]=0.0;
            spins::Sz[i]=0.0;
        }
        std::ifstream ifs("vampire.in");
        std::string dump;
        if(!ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Cannot read the vampire spin config file, check it is called vampire.in");
        }
        else
        {
            //get rid of the first line
            std::getline(ifs,dump);
        }
        for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
                {
                    for(unsigned int t = 0 ; t < geom::ucm.NumAtomsUnitCell() ; t++)
                    {
                        //get the atom number
                        int an=geom::coords(i*geom::Nk[0]+geom::ucm.GetCoord(t,0),j*geom::Nk[1]+geom::ucm.GetCoord(t,1),k*geom::Nk[2]+geom::ucm.GetCoord(t,2),0);
                        if(an<0)
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("There is an inconsistency between the input (vampire) unit cell/dimensions and your input.");
                        }
                        //vampire outputs each atom in the unit cell so you have to check
                        //that the order is

                        ifs >> dump >> dump >> dump >> dump >> spins::Sx[an] >> spins::Sy[an] >> spins::Sz[an];
//                        std::cout << spins::Sx[an] << "\t" << spins::Sy[an] << "\t" << spins::Sz[an] << std::endl;
                    }
                }
            }
        }
    }

}
