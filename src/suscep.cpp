// File: suscep.h
// Author: Tom Ostler
// Created: 25 Jan 2013
// Last-modified: 27 Mar 2013 17:26:49
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "../inc/arrays.h"
#include "../inc/defines.h"
#include "../inc/util.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/intmat.h"
#include "../inc/fields.h"
#include "../inc/llg.h"
#include "../inc/sim.h"
void sim::suscep(int argc,char *argv[])
{
	config::printline(config::Info);
	config::Info.width(45);config::Info << std::right << "*" << "**Chi(T) details***" << std::endl;
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

	libconfig::Setting &setting = config::cfg.lookup("suscep");
	double lT=0.0,uT=0.0,dT=0.0,convmean=0.0,convvar=0.0,met=0.0,et=0.0;
	setting.lookupValue("lower_temp",lT);
	FIXOUT(config::Info,"Lower temperature:" << lT << std::endl);
	setting.lookupValue("upper_temp",uT);
	FIXOUT(config::Info,"Upper temperature:" << uT << std::endl);
	setting.lookupValue("temp_step",dT);
	FIXOUT(config::Info,"Temperature step:" << dT << std::endl);
	setting.lookupValue("mean_tolerance",convmean);
	setting.lookupValue("variance_tolerance",convvar);
	FIXOUT(config::Info,"Converging mean to:" << convmean << std::endl);
	FIXOUT(config::Info,"Converging variance to:" << convvar << std::endl);
	setting.lookupValue("MaxRunTime",met);
	setting.lookupValue("EquilTime",et);
	FIXOUT(config::Info,"Maximum run time per temperature:" << met << " seconds" << std::endl);
	unsigned int mrts=int(met/llg::dt),ets=int(et/llg::dt);
	FIXOUT(config::Info,"Maximum run timesteps:" << mrts << std::endl);
	std::string opf;
	setting.lookupValue("SuscepFile",opf);
	FIXOUT(config::Info,"Outputting susceptibility data to:" << opf << std::endl);
    bool outmag=false;
    std::ofstream magout;
    std::string MagFilestr;
    setting.lookupValue("outmag",outmag);
    FIXOUT(config::Info,"Output magnetization data for each temperature:" << config::isTF(outmag) << std::endl);
    if(outmag==true)
    {
        setting.lookupValue("MagFileName",MagFilestr);
        FIXOUT(config::Info,"Magnetization files output to:" << MagFilestr << std::endl);
    }

	std::ofstream ofs(opf.c_str());
	if(!ofs.is_open())
	{
		error::errPreamble(__FILE__,__LINE__);
		error::errMessage("Could not open file stream for outputting magnetization data.");
	}
	else
	{
		ofs << "#Temperature\tMean" << std::endl;
	}
	util::RunningStat MS,mzr,mz2r,mxr,mx2r,myr,my2r;
	double mzro=0.0,mz2ro=0.0,mxro=0.0,mx2ro=0.0,myro=0.0,my2ro=0.0;
	//temperature loop
	for(double T = lT ; T < uT ; T+=dT)
	{
		config::printline(config::Info);
		FIXOUT(config::Info,"Converging temperature:" << T << std::endl);
		llg::T=T;
		MS.Clear();
		double oldmean=0.0;
		bool convTF=false;
		if(outmag)
		{
			std::stringstream MagFilesstr;
			MagFilesstr << MagFilestr << "_" << T << ".dat";
			std::string mopf=MagFilesstr.str();
			magout.open(mopf.c_str());
			if(!magout.is_open())
			{
				error::errPreamble(__FILE__,__LINE__);
				error::errMessage("Could not open file for outputting magnetization data");
			}
		}

		for(unsigned int t = 0 ; t < ets ; t++)
		{
			llg::integrate(t);
		}
		for(unsigned int t = ets ; t < mrts+ets ; t++)
		{
			llg::integrate(t);
			if(t%spins::update==0)
			{
				const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
				const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
				const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
				const double modm=sqrt(mx*mx+my*my+mz*mz);

				MS.Push(modm);
				//longitudinal components
				mzr.Push(mz);
				mz2r.Push(mz*mz);
				//transverse components
				mxr.Push(mx);
				myr.Push(my);
				mx2r.Push(mx*mx);
				my2r.Push(my*my);
				magout << T << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
				config::Info.width(15);config::Info << "| Mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS.Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MS.Mean()-oldmean) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS.Variance() << " [ " << convvar << " ]|" << std::endl;
				ofs.width(15);ofs << "| <m_z> = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mzr.Mean() << " | delta = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(mzr.Mean()-mzro) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mzr.Variance() << " [ " << convvar << " ]|" << std::endl;
				ofs.width(15);ofs << "| <m_z^2> = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mz2r.Mean() << " | delta = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(mz2r.Mean()-mz2ro) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mz2r.Variance() << " [ " << convvar << " ]|" << std::endl;
				ofs.width(15);ofs << "| <m_x> = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mxr.Mean() << " | delta = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(mxr.Mean()-mxro) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mxr.Variance() << " [ " << convvar << " ]|" << std::endl;
				ofs.width(15);ofs << "| <m_x^2> = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mx2r.Mean() << " | delta = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(mx2r.Mean()-mx2ro) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << mx2r.Variance() << " [ " << convvar << " ]|" << std::endl;
				ofs.width(15);ofs << "| <m_y> = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << myr.Mean() << " | delta = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(myr.Mean()-myro) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << myr.Variance() << " [ " << convvar << " ]|" << std::endl;
				ofs.width(15);ofs << "| <m_y^2> = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << my2r.Mean() << " | delta = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(my2r.Mean()-my2ro) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << my2r.Variance() << " [ " << convvar << " ]|" << std::endl;
				config::printline(ofs);






				if(((fabs(MS.Mean()-oldmean)) < convmean) && (MS.Variance()<convvar) && ((fabs(mzr.Mean()-mzro)) < convmean) && (mzr.Variance()<convvar) && ((fabs(mxr.Mean()-mxro)) < convmean) && (mxr.Variance()<convvar) && ((fabs(myr.Mean()-myro)) < convmean) && (myr.Variance()<convvar) && ((fabs(mz2r.Mean()-mz2ro)) < convmean) && (mz2r.Variance()<convvar) && ((fabs(mx2r.Mean()-mx2ro)) < convmean) && (mx2r.Variance()<convvar) && ((fabs(my2r.Mean()-my2ro)) < convmean) && (my2r.Variance()<convvar))
				{
					convTF=true;
					break;
				}
				oldmean=MS.Mean();
				mzro=mzr.Mean();
				mz2ro=mz2r.Mean();
				mxro=mxr.Mean();
				mx2ro=mx2r.Mean();
				myro=myr.Mean();
				my2ro=my2r.Mean();

			}
		}
		ofs << T << "\t" << MS.Mean() << "\t" << (mat::mu*geom::nspins/(2.0*1.38e-23*llg::T))*(mz2r.Mean()-mzr.Mean()*mzr.Mean()) << "\t" << (mat::mu*geom::nspins/(1.38e-23*llg::T))*(mx2r.Mean()+my2r.Mean()-mxr.Mean()*mxr.Mean()-myr.Mean()*myr.Mean()) << std::endl;

		FIXOUT(config::Info,"Converged?" << config::isTF(convTF) << std::endl);
	}
	ofs.close();
	if(ofs.is_open())
	{
		error::errPreamble(__FILE__,__LINE__);
		error::errWarning("Could not close file for outputting magnetization data.");
	}
}
