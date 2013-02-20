// File: mvt.h
// Author: Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 20 Feb 2013 13:01:27
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
    std::string order_param;
	double lT=0.0,uT=0.0,dT=0.0,convmean=0.0,convvar=0.0,met=0.0,minrt=0.0;
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
    setting.lookupValue("MinRunTime",minrt);
    FIXOUT(config::Info,"Minimum runtime per temp step:" << minrt << std::endl);
	unsigned int mrts=int(met/llg::dt),eminrts=int(minrt/llg::dt);
	FIXOUT(config::Info,"Maximum run timesteps:" << mrts << std::endl);
    FIXOUT(config::Info,"Min rum timesteps:" << eminrts << std::endl);
	std::string opf;
	setting.lookupValue("MvTFile",opf);
	FIXOUT(config::Info,"Outputting magnetization data to:" << opf << std::endl);
    setting.lookupValue("order_param",order_param);
    FIXOUT(config::Info,"Outputting order parameter:" << order_param << std::endl);
    unsigned int nslat=1;
    std::string afmtype;
    double *mx=NULL,*my=NULL,*mz=NULL,*modm=NULL,*ncount=NULL;;
    unsigned int *latlookup=NULL;
    setting.lookupValue("Number_sublattices",nslat);
    FIXOUT(config::Info,"Number of sublattices in order parameter calculation:" << nslat << std::endl);
    try
    {
        mx = new double[nslat];
        my = new double[nslat];
        mz = new double[nslat];
        modm = new double[nslat];
        ncount = new double[nslat];
        latlookup = new unsigned int[geom::nspins];
        for(unsigned int i = 0 ; i < nslat ; i++)
        {
            mx[i]=0;
            my[i]=0;
            mz[i]=0;
            modm[i]=0;
            ncount[i]=0;
        }
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int zc=geom::lu(i,2);
            latlookup[i]=zc%nslat;
            ncount[zc%nslat]++;
        }
        unsigned int tc=0;
        for(unsigned int i = 0 ; i < nslat ; i++)
        {
            tc+=ncount[i];
        }
        if(tc!=geom::nspins)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Mismatch in counting sublattices spins");
        }
    }
    catch(...)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not allocate staggered magnetization memory");
    }
    if(nslat>1)
    {
        setting.lookupValue("AFM_type",afmtype);
        FIXOUT(config::Info,"Type of AFM expected:" << afmtype << std::endl);
        if(afmtype!="layered")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Type of AFM not currently coded");
        }
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
	util::RunningStat MS[nslat];
	util::RunningStat MSx[nslat],MSy[nslat],MSz[nslat];
    bool conv[nslat];
	//temperature loop
	for(double T = lT ; T < uT ; T+=dT)
	{
		config::printline(config::Info);
		FIXOUT(config::Info,"Converging temperature:" << T << std::endl);
		llg::T=T;
		double oldmean[nslat];
        for(unsigned int nl = 0 ; nl < nslat ; nl++)
        {
            oldmean[nl]=0.0;
            conv[nl]=false;
            MS[nl].Clear();
			MSx[nl].Clear();
			MSy[nl].Clear();
			MSz[nl].Clear();
        }
		bool convTF=false;
		for(unsigned int t = 0 ; t < eminrts ; t++)
		{
			llg::integrate(t);
		}
		for(unsigned int t = 0 ; t < mrts ; t++)
		{
			llg::integrate(t);
			if(t%spins::update==0)
			{
                for(unsigned int i = 0 ; i < nslat ; i++)
                {
                    mx[i]=0.0;
                    my[i]=0.0;
                    mz[i]=0.0;
                    modm[i]=0.0;
                }

                if(nslat==1)
                {
                    mx[0] = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
                    my[0] = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
                    mz[0] = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
                    modm[0]=sqrt(mx[0]*mx[0]+my[0]*my[0]+mz[0]*mz[0]);
                }
                else
                {
                    for(unsigned int i = 0 ; i < geom::nspins ; i++)
                    {
                        //sublattice owner
                        unsigned int slowner=latlookup[i];
                        mx[slowner]+=spins::Sx[i];
                        my[slowner]+=spins::Sy[i];
                        mz[slowner]+=spins::Sz[i];
                    }
                    for(unsigned int i = 0 ; i < nslat ; i++)
                    {
                        mx[i]/=ncount[i];
                        my[i]/=ncount[i];
                        mz[i]/=ncount[i];
                        modm[i]=sqrt(mx[i]*mx[i]+my[i]*my[i]+mz[i]*mz[i]);
                    }
                }





				if(t>int(10e-12/llg::dt))
				{
                    for(unsigned int i = 0 ; i < nslat ; i++)
                    {
                        MS[i].Push(modm[i]);
						MSx[i].Push(mx[i]);
						MSy[i].Push(my[i]);
						MSz[i].Push(mz[i]);
                        config::Info.width(15);config::Info << "Sublattice " << i << "| Mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS[i].Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MS[i].Mean()-oldmean[i]) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS[i].Variance() << " [ " << convvar << " ]|" << std::endl;
                        if(((fabs(MS[i].Mean()-oldmean[i])) < convmean) && (MS[i].Variance()<convvar) && t > int(75e-12/llg::dt))
                        {
                            conv[i]=true;
                        }
                        oldmean[i]=MS[i].Mean();
                    }
                    convTF=true;
                    for(unsigned int i = 0 ; i < nslat ; i++)
                    {
                        if(conv[i]==false)
                        {
                            convTF=false;
                        }
                    }
                    if(convTF==true)
                    {
                        break;
                    }

				}

			}
/*			std::cout << t;
			for(unsigned int i = 0 ; i < nslat ; i++)
			{
				std::cout << "\t" << mx[i] << "\t" << my[i] << "\t" << mz[i] << "\t" << modm[i];
			}
			std::cout << std::endl;*/
		}
        ofs << T;
        for(unsigned int i = 0 ; i < nslat ; i++)
        {
            ofs << "\t" << MSx[i].Mean() << "\t" << MSy[i].Mean() << "\t" << MSz[i].Mean() << "\t" << MS[i].Mean();
        }
		ofs << std::endl;


		FIXOUT(config::Info,"Converged?" << config::isTF(convTF) << std::endl);
	}
	ofs.close();
	if(ofs.is_open())
	{
		error::errPreamble(__FILE__,__LINE__);
		error::errWarning("Could not close file for outputting magnetization data.");
	}

	delete [] mx;
	delete [] my;
	delete [] mz;
	delete [] modm;
	delete [] ncount;
	delete [] latlookup;
	mx=NULL;my=NULL;mz=NULL;modm=NULL;ncount=NULL;latlookup=NULL;


}
