// File: mvt.h
// Author: Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 28 Jan 2013 10:42:59
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
void sim::MvT(int argc,char *argv[])
{
	config::printline(config::Info);
	config::Info.width(45);config::Info << std::right << "*" << "**M(T) details***" << std::endl;
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

	libconfig::Setting &setting = config::cfg.lookup("mvt");
	double lT=0.0,uT=0.0,dT=0.0,convmean=0.0,convvar=0.0,met=0.0;
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
	FIXOUT(config::Info,"Maximum run time per temperature:" << met << " seconds" << std::endl);
	unsigned int mrts=int(met/llg::dt);
	FIXOUT(config::Info,"Maximum run timesteps:" << mrts << std::endl);
	std::string opf;
	setting.lookupValue("MvTFile",opf);
	FIXOUT(config::Info,"Outputting magnetization data to:" << opf << std::endl);
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
	util::RunningStat MS;
	//temperature loop
	for(double T = lT ; T < uT ; T+=dT)
	{
        double modm=0.0;
		config::printline(config::Info);
		FIXOUT(config::Info,"Converging temperature:" << T << std::endl);
		llg::T=T;
		MS.Clear();
		double oldmean=0.0;
		bool convTF=false;
		for(unsigned int t = 0 ; t < mrts ; t++)
		{
			llg::integrate(t);
			if(t%spins::update==0)
			{
				const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
				const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
				const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
				modm=sqrt(mx*mx+my*my+mz*mz);
				if(t>int(50e-12/llg::dt))
				{
					MS.Push(modm);
					config::Info.width(15);config::Info << "| Mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS.Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MS.Mean()-oldmean) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS.Variance() << " [ " << convvar << " ]|" << std::endl;
					if(((fabs(MS.Mean()-oldmean)) < convmean) && (MS.Variance()<convvar) && t > int(75e-12/llg::dt))
					{
						convTF=true;
						break;
					}
					oldmean=MS.Mean();
				}

			}
		}
        ofs << T << "\t" << modm << std::endl;


		FIXOUT(config::Info,"Converged?" << config::isTF(convTF) << std::endl);
	}
	ofs.close();
	if(ofs.is_open())
	{
		error::errPreamble(__FILE__,__LINE__);
		error::errWarning("Could not close file for outputting magnetization data.");
	}




}
