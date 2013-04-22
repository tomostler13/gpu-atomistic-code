// File: llg.cpp
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 22 Apr 2013 14:40:14
#include "../inc/llg.h"
#include "../inc/llgCPU.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/mat.h"
#include "../inc/spins.h"
#include "../inc/llg.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/exch.h"
#include <cmath>
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace llg
{
    double applied[3]={0,0,0},T,dt,rdt;
    Array<double> osT,llgpf;
    //real space correlation function
    bool rscf=false;
    //static structure factor
    bool ssf=false;
    //on site applied field?
    bool osHapp=false;
    //on site temperature?
    bool osTemp=false;
    //type of on site applied field
    std::string osk;

    std::string rscfstr;
    std::ofstream rscfs;
    std::string ssffs;
    std::ofstream ssffstr;
	void initLLG(int argc,char *argv[])
	{
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**LLG details***" << std::endl;

        SUCCESS(config::Info);
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
        llgpf.resize(geom::nspins);
        libconfig::Setting &setting = config::cfg.lookup("llg");
        setting.lookupValue("dt",dt);
        setting.lookupValue("RealSpaceCorrelations",rscf);
        setting.lookupValue("StaticStructureFactor",ssf);

        setting.lookupValue("RSCFile",rscfstr);
        setting.lookupValue("SSFile",ssffs);
        setting.lookupValue("Onsite_Applied",osHapp);
        FIXOUT(config::Info,"Onsite applied field?:" << config::isTF(osHapp) << std::endl);
        FIXOUT(config::Info,"Resizing applied field arrays:" << std::flush);
        setting.lookupValue("Onsite_Temp",osTemp);
        if(osTemp)
        {
            osT.resize(geom::nspins);
            osT.IFill(0);
        }
        FIXOUT(config::Info,"Onsite temperature?:" << config::isTF(osTemp) << std::endl);
        if(osHapp)
        {
            fields::HAppx.resize(geom::nspins);
            fields::HAppy.resize(geom::nspins);
            fields::HAppz.resize(geom::nspins);
        }
        SUCCESS(config::Info);
        setting.lookupValue("Onsite_Kind",osk);
        FIXOUT(config::Info,"Kind of onsite applied field:" << osk);
        FIXOUT(config::Info,"Outputting correlation functions to file:" << rscfstr << std::endl);
        FIXOUT(config::Info,"Calculating real space correlation functions:" << config::isTF(rscf) << std::endl);
        if(ssf)
        {
            FIXOUT(config::Info,"Calculating static structure factor:" << config::isTF(ssf) << std::endl);
            FIXOUT(config::Info,"Outputting to file:" << ssffs << std::endl);
        }
        for(unsigned int i = 0 ; i < 3 ;i++)
        {
            applied[i]=setting["applied"][i];
        }
        FIXOUTVEC(config::Info,"Applied field:",applied[0],applied[1],applied[2]);
        FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
        rdt=dt*mat::gamma;
		FIXOUT(config::Info,"Reduced timestep:" << rdt << std::endl);
        FIXOUT(config::Info,"Setting initial sigma" << std::flush);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            mat::sigma[i] = sqrt(2.0*1.38e-23*mat::lambda[i]/(mat::mu[i]*mat::muB*dt*mat::gamma));
            llgpf[i] = -1./(1.0+mat::lambda[i]*mat::lambda[i]);
        }
        SUCCESS(config::Info);
//        FIXOUT(config::Info,"Sigma prefactor:" << mat::sigma << std::endl);
//		FIXOUT(config::Info,"Prefactor to LLG equation:" << llgpf << std::endl);

	}
    void integrate(unsigned int& t)
    {
        #ifdef CUDA
        if(config::useintmat)
        {
            if(osTemp)
            {
                cullg::llgGPU(t,osT);
            }
            else
            {
                cullg::llgGPU(t);
            }
        }
        else
        {
            if(osTemp)
            {
                cullg::llgGPU(t,osT,exch::xadj,exch::adjncy);
            }
            else
            {
                cullg::llgGPU(t,exch::xadj,exch::adjncy);
            }
        }
        #else
        if(config::useintmat)
        {
            if(osTemp)
            {
                llgCPU::llgCPU(t,osT);
            }
            else
            {
                llgCPU::llgCPU(t);
            }
        }
        else
        {
            if(osTemp)
            {
                llgCPU::llgCPU(t,osT,exch::xadj,exch::adjncy);
            }
            else
            {
                llgCPU::llgCPU(t,exch::xadj,exch::adjncy);
            }
        }
        #endif
		if(t%spins::update==0)
		{
            if(rscf)
            {
                spins::calcRealSpaceCorrelationFunction(t);
            }
            if(ssf)
            {
                spins::calcStaticStructureFactor(t);
            }
		}

    }

}
