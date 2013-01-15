// File: mat.cpp
// Author:Tom Ostler
// Last-modified: 19 Dec 2012 10:12:55
#include "../inc/mat.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include <libconfig.h++>
#include <cassert>
#include <sstream>
#include <string>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace mat
{
    //initialization checker
    bool mi=false;
	//spontaneous magnetization (J/T/m3)
	double Ms=0.0;
	//bath coupling constant (no units)
	double lambda=0.0;
	//Curie temperature (K)
	double Tc=0.0;
	//gyromagnetic ratio
	double gamma=0.0;
    //Anisotropy constant
    double K=0.0;
    //a prefactor for calculating perpendicular anisotropy
    //from MFA
    double chippf=0.0;
    //Number of sublattices
    unsigned int nlat=0;
    //Localized atomic magnetic moment of sublattice i
    Array<double> mu;
    //concentrations of the sublattices
    Array<double> conc;
    //string to identify the material, e.g. FePt, Fe, GdFeCo etc.
    std::string materialid;
	void initMat(int argc,char *argv[])
	{

		assert(geom::gi);
		config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Material details***" << std::endl;
		config::cfg.readFile(argv[1]);
		libconfig::Setting &setting = config::cfg.lookup("material");
		try
		{
			setting.lookupValue("Ms",Ms);
			setting.lookupValue("lambda",lambda);
			setting.lookupValue("gamma",gamma);
			setting.lookupValue("Tc",Tc);
            setting.lookupValue("K",K);
            setting.lookupValue("nlat",nlat);
            mu.resize(nlat);
            conc.resize(nlat);
            if(nlat<1)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Number of sublattices not recognized");
            }
            double cc=0;
            for(unsigned int i = 0 ; i < nlat ; i++)
            {
                mu[i]=setting["mu"][i];
                conc[i]=setting["conc"][i];
                cc+=conc[i];
            }
            if(cc<(1-1e-12) && cc>(1.0+1e-12))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Concentrations of sublattices not set correctly");
            }
            setting.lookupValue("materialid",materialid);
		}
		catch(const libconfig::SettingTypeException &stex)
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Setting type error");
		}
        FIXOUT(config::Info,"Material identified as:" << materialid << std::endl);
		FIXOUT(config::Info,"Ms:" << Ms << " J/T/m^3" << std::endl);
		FIXOUT(config::Info,"lambda:" << lambda << std::endl);
		FIXOUT(config::Info,"gamma:" << gamma << "J^{-1}T^{-1}" << std::endl);
		FIXOUT(config::Info,"Curie temperature:" << Tc << " K" << std::endl);
        FIXOUT(config::Info,"Anisotropy constant (J/m^3)" << K << std::endl);
        FIXOUT(config::Info,"Number of sublattices:" << nlat << std::endl);

        chippf=Ms/(2.0*K);
        for(unsigned int i = 0 ; i < nlat ; i++)
        {
            std::stringstream sstr;
            sstr << "Concentration of sublattice " << i << ":";
            std::string str=sstr.str();
            config::Info.width(75);config::Info << std::left << str << conc[i] << std::endl;
        }

        for(unsigned int i = 0 ; i < nlat ; i++)
        {
            std::stringstream sstr;
            sstr << "Magnetic moment of sublattice " << i << ":";
            std::string str=sstr.str();

            config::Info.width(75);config::Info << std::left << str << conc[i] << std::endl;
        }

		assert(Ms>0);
		assert(Tc>0);
		assert(gamma>0);
        mi=true;
	}

}
