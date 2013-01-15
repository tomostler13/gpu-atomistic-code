// File: fields.cpp
// Author:Tom Ostler
// Last-modified: 13 Dec 2012 19:47:40
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../inc/array3d.h"
#include "../inc/array.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/tdp.h"
#include "../inc/mf.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace fields
{

	Array3D<fftw_complex> Hkx;
	Array3D<fftw_complex> Hky;
	Array3D<fftw_complex> Hkz;

    Array<double> Hx;
    Array<double> Hy;
    Array<double> Hz;

    Array<double> GW1x;
    Array<double> GW1y;
    Array<double> GW1z;

    Array<double> GW2x;
    Array<double> GW2y;
    Array<double> GW2z;

    bool checkthermal=false;
    bool checkdipolar=false;
    bool checkexchange=false;

    //how often to update the dipolar field
    unsigned int dfu=1;
    unsigned int anist=0;
    //function pointer to the correct form of anisotropy
    //allows the anisotropy form to be determined at runtime
    void (*anisfp)(double *,const double*,unsigned int);

    //applied field details
    unsigned int ff;
    //function pointer to applied field form
    void (*appliedfp)(double *);
    double fH[3]={0,0,0};

	std::string dfc;
	bool fi=false;
	fftw_plan HxP=NULL,HyP=NULL,HzP=NULL;

	void initFields(int argc,char *argv[])
	{
		//make sure we have specified the system
		assert(geom::gi);
        assert(mf::mfi);
        GW1x.resize(geom::ss);
        GW1y.resize(geom::ss);
        GW1z.resize(geom::ss);
        GW2x.resize(geom::ss);
        GW2y.resize(geom::ss);
        GW2z.resize(geom::ss);
        GW1x.IFill(0);
        GW1y.IFill(0);
        GW1z.IFill(0);
        GW2x.IFill(0);
        GW2y.IFill(0);
        GW2z.IFill(0);
		config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Field details***" << std::endl;
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

		libconfig::Setting &setting = config::cfg.lookup("fields");
		try
		{
			setting.lookupValue("dipfieldcalc",dfc);
            setting.lookupValue("dipupdate",dfu);
            setting.lookupValue("anisotropy",anist);
            setting.lookupValue("checkthermal",checkthermal);
            setting.lookupValue("checkdipolar",checkdipolar);
            setting.lookupValue("checkexchange",checkexchange);
            if(anist==0)
            {
                anisfp = anis0;
            }
            else if(anist==1)
            {
                anisfp = anis1u;
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Anisotropy type not recognized");
            }
            setting.lookupValue("fieldform",ff);
            FIXOUT(config::Info,"Thermal fields:" << config::isTF(checkthermal) << std::endl);
            FIXOUT(config::Info,"Dipolar fields:" << config::isTF(checkdipolar) << std::endl);
            FIXOUT(config::Info,"Exchange fields:" << config::isTF(checkexchange) << std::endl);


            if(ff==0)
            {

                FIXOUT(config::Info,"Applied field type:" << "Fixed" << std::endl);
                for(unsigned int i = 0 ; i < 3 ; i++)
                {
                    fH[i]=setting["appliedfixed"][i];
                }
                appliedfp = applied0;
                FIXOUT(config::Info,"Field:" << "(" << fH[0] << "," << fH[1] << "," << fH[2] << ")" << std::endl);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Applied field not recognised");
            }
		}
		catch(const libconfig::SettingTypeException &stex)
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Setting type error");
		}
        std::string str;
        if(anist==0)
        {
            str = "No anisotropy";
        }
        else if(anist==1)
        {
            str = "-(m_x*e_x + m_y*e_y)/chiperp";
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Anisotropy type not recognized");
        }

        FIXOUT(config::Info,"Anisotropy type:" << str << std::endl);
		FIXOUT(config::Info,"Method for calculating dipolar fields:" << dfc << std::endl);
		FIXOUT(config::Info,"Resizing field arrays:" << std::flush);
		Hkx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Hky.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
		Hkz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);

        Hx.resize(geom::ss);
        Hy.resize(geom::ss);
        Hz.resize(geom::ss);


		Hkx.IFill(0);
        Hky.IFill(0);
        Hkz.IFill(0);
        Hx.IFill(0);
        Hy.IFill(0);
        Hz.IFill(0);
		config::Info << "Done" << std::endl;
		FIXOUT(config::Info,"Planning backward transform of fields:" << std::flush);
        //plan the transforms as in-place as we do not need to use the fft arrays
        //as we copy the data back to the normal field arrayl
		HxP = fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Hkx.ptr(),Hkx.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
		HyP = fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Hky.ptr(),Hky.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
		HzP = fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Hkz.ptr(),Hkz.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);


        atexit(destroyFieldPlans);
		config::Info << "Done" << std::endl;
		fi=true;
	}

    void transformFieldsBack()
    {
        assert(HxP!=NULL);
        assert(HyP!=NULL);
        assert(HzP!=NULL);
        assert(fi);
        fftw_execute(HxP);
        fftw_execute(HyP);
        fftw_execute(HzP);
        unsigned int count=0;
        //normalize the fourier transform
        for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                {
                    Hkx(i,j,k)[0]/=double(geom::zps);
                    Hky(i,j,k)[0]/=double(geom::zps);
                    Hkz(i,j,k)[0]/=double(geom::zps);
                }
            }
        }
        //copy fields to normal arrays
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            unsigned int ijk[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Hx(i)=Hkx(ijk[0],ijk[1],ijk[2])[0];
            Hy(i)=Hky(ijk[0],ijk[1],ijk[2])[0];
            Hz(i)=Hkz(ijk[0],ijk[1],ijk[2])[0];
        }


    }

    void destroyFieldPlans(void)
    {
        fftw_destroy_plan(HxP);
        fftw_destroy_plan(HyP);
        fftw_destroy_plan(HzP);
        HxP=NULL;
        HyP=NULL;
        HzP=NULL;
        assert(HxP==NULL);
        assert(HyP==NULL);
        assert(HzP==NULL);
    }

    void bfdip()
    {
        assert(spins::si);
        fields::Hx.IFill(0);
        fields::Hy.IFill(0);
        fields::Hz.IFill(0);
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            double h[3]={0,0,0};
            for(unsigned int j = 0 ; j < geom::ss ; j++)
            {
                if(i!=j)
                {
                    double rij[3]={(double(geom::lu(j,0))-double(geom::lu(i,0)))*geom::gs[0],(double(geom::lu(j,1))-double(geom::lu(i,1)))*geom::gs[1],(double(geom::lu(j,2))-double(geom::lu(i,2)))*geom::gs[2]};
                    double mrij=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                    double oomrij3=1./(mrij*mrij*mrij);
                    double eij[3]={rij[0]/mrij,rij[1]/mrij,rij[2]/mrij};
                    double sj[3]={spins::sx(j),spins::sy(j),spins::sz(j)};
                    double sjdote=sj[0]*eij[0]+sj[1]*eij[1]+sj[2]*eij[2];
                    h[0]+=(3.0*sjdote*eij[0]-sj[0])*oomrij3;
                    h[1]+=(3.0*sjdote*eij[1]-sj[1])*oomrij3;
                    h[2]+=(3.0*sjdote*eij[2]-sj[2])*oomrij3;
                }
            }
            fields::Hx[i]=h[0]*1e-7*geom::gsV*mat::Ms;
            fields::Hy[i]=h[1]*1e-7*geom::gsV*mat::Ms;
            fields::Hz[i]=h[2]*1e-7*geom::gsV*mat::Ms;
        }
    }

    void anis0(double *h,const double *s,unsigned int i)
    {
        h[0]+=0;
        h[1]+=0;
        h[2]+=0;
    }
    void anis1u(double *h,const double *s,unsigned int i)
    {
        h[0]+=-s[0]/tdp::syschiperp[i];
        h[1]+=-s[1]/tdp::syschiperp[i];
    }
    void applied0(double *h)
    {
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            h[i]+=fH[i];
        }
    }
}
