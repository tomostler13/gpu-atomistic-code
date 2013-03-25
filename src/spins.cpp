// File: spins.cpp
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 25 Mar 2013 17:50:57
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/mat.h"
#include "../inc/intmat.h"
#include "../inc/defines.h"
#include "../inc/maths.h"
#include "../inc/llg.h"
namespace spins
{
    Array3D<fftw_complex> Skx,Sky,Skz;
    Array3D<double> Srx,Sry,Srz;
    Array3D<double> Sznzp;
    Array3D<fftw_complex> Sqznzp;
    unsigned int nzpcplxdim=0;
    double normsize=0;

    Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
    fftw_plan SxP,SyP,SzP,SzcfPF,SzcfPB;
	unsigned int update=0;
	std::ifstream sfs;
    void initSpins(int argc,char *argv[])
    {
		config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Spin details***" << std::endl;
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
		FIXOUT(config::Info,"Resizing arrays:" << std::flush);
        Skx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Sky.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Skz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Srx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Sry.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Srz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Srx.IFill(0);
        Sry.IFill(0);
        Srz.IFill(0);

        Sx.resize(geom::nspins);
        Sy.resize(geom::nspins);
        Sz.resize(geom::nspins);
		SUCCESS(config::Info);
		std::string sc;
		libconfig::Setting &setting = config::cfg.lookup("spins");
        if(llg::rscf)
        {
            Sznzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);

            nzpcplxdim=(geom::dim[2]*geom::Nk[2]/2)+1;
            Sqznzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
            SzcfPF = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sznzp.ptr(),Sqznzp.ptr(),FFTW_ESTIMATE);
            SzcfPB = fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sqznzp.ptr(),Sznzp.ptr(),FFTW_ESTIMATE);
            normsize=geom::dim[0]*geom::Nk[0]*geom::dim[0]*geom::Nk[0];
            normsize*=geom::dim[1]*geom::Nk[1]*geom::dim[1]*geom::Nk[1];
            normsize*=geom::dim[2]*geom::Nk[2]*geom::dim[2]*geom::Nk[2];
        }


		setting.lookupValue("update",update);
		FIXOUT(config::Info,"Spin update:" << update << " (timesteps)" << std::endl);
		if(update<1)
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Spin data should be updated on cpu, i.e. spins::update>0");
		}
		setting.lookupValue("spinconfig",sc);
		FIXOUT(config::Info,"Initial spin config specifier method:" << sc << std::endl);
		if(sc=="file")
		{
			std::string spinfile;
			try
			{
				setting.lookupValue("sf",spinfile);
			}
			catch(const libconfig::SettingTypeException &stex)
			{
				error::errPreamble(__FILE__,__LINE__);
				error::errMessage("Setting type error");
			}
			FIXOUT(config::Info,"Reading spin data from file:" << spinfile << std::flush);
			sfs.open(spinfile.c_str());
			if(!sfs.is_open())
			{
				error::errPreamble(__FILE__,__LINE__);
				error::errMessage("Could not open file for reading spin data");
			}
			else
			{
				for(unsigned int i = 0 ; i < geom::nspins ; i++)
				{
                    sfs >> Sx(i);
                    sfs >> Sy(i);
                    sfs >> Sz(i);
				}
				sfs.close();
				if(sfs.is_open())
				{
					error::errPreamble(__FILE__,__LINE__);
					error::errWarning("Could not close spin file for reading in spin data");
				}
			}
		}
		else if(sc=="align")
		{
			double sr[3]={0,0,0};
			try
			{
				for(unsigned int l = 0 ; l < 3 ; l++){sr[l]=setting["sc"][l];}
			}
			catch(const libconfig::SettingTypeException &stex)
			{
				error::errPreamble(__FILE__,__LINE__);
				error::errMessage("Setting type error");
			}
			FIXOUT(config::Info,"Spins aligned to:" << "(" << sr[0] << "," << sr[1] << "," << sr[2] << ")" <<std::endl);
			for(unsigned int i = 0 ; i < geom::nspins ; i++)
			{
                Sx(i)=sr[0];
                Sy(i)=sr[1];
                Sz(i)=sr[2];//sqrt(1.0-(Sx[i]*Sx[i]+Sy[i]*Sy[i]));//+sr[2];
			}
		}
		else if(sc=="random")
		{
			for(unsigned int i = 0 ; i < geom::nspins ; i++)
			{
                maths::sphereRandom(Sx(i),Sy(i),Sz(i));
			}
		}
		else
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Spin configuration not recognised, check you config file.");
		}

		FIXOUT(config::Info,"Planning r2c and c2r transforms:" << std::flush);
        //forward transform of spin arrays
        SxP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Srx.ptr(),Skx.ptr(),FFTW_ESTIMATE);
        SyP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Sry.ptr(),Sky.ptr(),FFTW_ESTIMATE);
        SzP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Srz.ptr(),Skz.ptr(),FFTW_ESTIMATE);
		SUCCESS(config::Info);

    }
    void FFTForward()
    {
        Skx.IFill(0);
        Sky.IFill(0);
        Skz.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Srx(lc[0],lc[1],lc[2])=Sx[i];
            Sry(lc[0],lc[1],lc[2])=Sy[i];
            Srz(lc[0],lc[1],lc[2])=Sz[i];
//            std::cout << "Lookup\t" << lc[0] << "\t" << lc[1] << "\t" << lc[2] << "\t" << Skx(lc[0],lc[1],lc[2])[0] << "\t" << Sky(lc[0],lc[1],lc[2])[0] << "\t" << Skz(lc[0],lc[1],lc[2])[0] << std::endl;
        }
        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
    }
    void eFFTForward()
    {
        Skx.IFill(0);
        Sky.IFill(0);
        Skz.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Srx(lc[0],lc[1],lc[2])=eSx[i];
            Sry(lc[0],lc[1],lc[2])=eSy[i];
            Srz(lc[0],lc[1],lc[2])=eSz[i];
//            std::cout << "Lookup\t" << lc[0] << "\t" << lc[1] << "\t" << lc[2] << "\t" << Skx(lc[0],lc[1],lc[2])[0] << "\t" << Sky(lc[0],lc[1],lc[2])[0] << "\t" << Skz(lc[0],lc[1],lc[2])[0] << std::endl;
        }
        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
    }
    void calcRealSpaceCorrelationFunction(unsigned int t)
    {
        Sqznzp.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Sznzp(lc[0],lc[1],lc[2])=Sz[i];
        }
        fftw_execute(SzcfPF);
        for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                {
                    double Sq[2]={Sqznzp(i,j,k)[0],Sqznzp(i,j,k)[1]};
                    double CCSq[2]={Sqznzp(i,j,k)[0],-Sqznzp(i,j,k)[1]};
                    Sqznzp(i,j,k)[0]=Sq[0]*CCSq[0]-Sq[1]*CCSq[1];//Sqznzp(i,j,k)[0]*Sqznzp(i,j,k)[0]-Sqznzp(i,j,k)[1]*Sqznzp(i,j,k)[1];
                    Sqznzp(i,j,k)[1]=Sq[0]*CCSq[1]+Sq[1]*CCSq[0];//Sqznzp(i,j,k)[0]*Sqznzp(i,j,k)[1]+Sqznzp(i,j,k)[1]*Sqznzp
                }
            }
        }
        fftw_execute(SzcfPB);
        std::stringstream sstr;
        sstr << llg::rscfstr << t;
        std::string tempstr=sstr.str();
        llg::rscfs.open(tempstr.c_str());
        if(!llg::rscfs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing correlation function");
        }
        for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i+=geom::Nk[0])
        {
            for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j+=geom::Nk[1])
            {
                //for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k+=geom::Nk[2])
                //{
                    int ijk[3]={i,j,0};
                    for(unsigned int c = 0 ; c < 3 ; c++)
                    {
                        if(ijk[c]>geom::dim[c]*geom::Nk[c]/2)
                        {
                            ijk[c]=ijk[c]-geom::dim[c]*geom::Nk[c];
                        }
                    }
                    llg::rscfs << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << "\t" << Sznzp(i,j,0)/normsize << std::endl;
                //}
            }
            llg::rscfs << std::endl;
        }
        llg::rscfs.close();
        if(llg::rscfs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close file for writing correlation function");
        }

    }

}
