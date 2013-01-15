// File: spins.cpp
// Author:Tom Ostler
// Last-modified: 04 Jan 2013 12:02:54
#include "../inc/array3d.h"
#include "../inc/spins.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/maths.h"
#include "../inc/geom.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace spins
{
	Array3D<fftw_complex> rsx;
	Array3D<fftw_complex> csx;
	Array3D<fftw_complex> rsy;
	Array3D<fftw_complex> csy;
	Array3D<fftw_complex> rsz;
	Array3D<fftw_complex> csz;
    Array<double> sx;
    Array<double> sy;
    Array<double> sz;
    Array<double> esx;
    Array<double> esy;
    Array<double> esz;

    //overall magnetization (reduced)
    double mx=0;
    double my=0;
    double mz=0;
    double modm=0;

	std::string sc;
	std::ifstream sfs;
	fftw_plan SxP=NULL;
    fftw_plan SyP=NULL;
    fftw_plan SzP=NULL;
	bool si=false;

	void initSpins(int argc,char *argv[])
	{
		assert(geom::gi);
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
		//resize the spin array
		rsx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        rsy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
		rsz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);

		csx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
		csy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
		csz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);

        sx.resize(geom::ss);
        sy.resize(geom::ss);
        sz.resize(geom::ss);

        esx.resize(geom::ss);
        esy.resize(geom::ss);
        esz.resize(geom::ss);


        sx.IFill(0);
        sy.IFill(0);
        sz.IFill(0);
        esx.IFill(0);
        esy.IFill(0);
        esz.IFill(0);

        csx.IFill(0);
        csy.IFill(0);
        csz.IFill(0);
        rsx.IFill(0);
        rsy.IFill(0);
        rsz.IFill(0);



		libconfig::Setting &setting = config::cfg.lookup("spins");
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
				for(unsigned int i = 0 ; i < geom::ss ; i++)
				{
                    sfs >> sx(i);
                    sfs >> sy(i);
                    sfs >> sz(i);
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
			for(unsigned int i = 0 ; i < geom::ss ; i++)
			{
                sx(i)=sr[0];
                sy(i)=sr[1];
                sz(i)=sr[2];
			}
		}
		else if(sc=="random")
		{
			for(unsigned int i = 0 ; i < geom::ss ; i++)
			{
                maths::sphereRandom(sx(i),sy(i),sz(i));
			}
		}
		else
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Spin configuration not recognised, check you config file.");
		}
		FIXOUT(config::Info,"Planning CPU (back and forward) transforms of spin array:" << std::flush);

        //forward transform of spin arrays
        SxP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],rsx.ptr(),csx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        SyP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],rsy.ptr(),csy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        SzP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],rsz.ptr(),csz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        atexit(destroySpinPlans);

		config::Info << "Done" << std::endl;
		si=true;
	}
    void destroySpinPlans(void)
    {
        fftw_destroy_plan(SxP);
        fftw_destroy_plan(SyP);
        fftw_destroy_plan(SzP);
        SxP=NULL;
        SyP=NULL;
        SzP=NULL;
        assert(SxP==NULL);
        assert(SyP==NULL);
        assert(SzP==NULL);
    }
    void forwardTransformSpinArrays()
    {
        assert(SxP!=NULL);
        assert(SyP!=NULL);
        assert(SzP!=NULL);
        csx.IFill(0);
        csy.IFill(0);
        csz.IFill(0);
        rsx.IFill(0);
        rsy.IFill(0);
        rsz.IFill(0);

        //copy the spin array to the fftw_complex array
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            unsigned int ijk[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            rsx(ijk[0],ijk[1],ijk[2])[0]=sx(i);
            rsy(ijk[0],ijk[1],ijk[2])[0]=sy(i);
            rsz(ijk[0],ijk[1],ijk[2])[0]=sz(i);
        }

        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
    }

}
