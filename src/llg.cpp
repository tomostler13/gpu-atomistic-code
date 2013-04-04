// File: llg.cpp
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 04 Apr 2013 20:12:16
#include "../inc/llg.h"
#include "../inc/llgCPU.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/mat.h"
#include "../inc/spins.h"
#include "../inc/llg.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include <cmath>
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace llg
{
    double applied[3]={0,0,0},T,dt,rdt,llgpf;
    //real space correlation function
    bool rscf=false;
    //on site applied field?
    bool osHapp=false;
    //type of on site applied field
    std::string osk;

    std::string rscfstr;
    std::ofstream rscfs;
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

        libconfig::Setting &setting = config::cfg.lookup("llg");
        setting.lookupValue("dt",dt);
        setting.lookupValue("RealSpaceCorrelations",rscf);

        setting.lookupValue("RSCFile",rscfstr);
        setting.lookupValue("Onsite_Applied",osHapp);
        FIXOUT(config::Info,"Onsite applied field?:" << config::isTF(osHapp) << std::endl);
        FIXOUT(config::Info,"Resizing applied field arrays:" << std::flush);
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
        for(unsigned int i = 0 ; i < 3 ;i++)
        {
            applied[i]=setting["applied"][i];
        }
        FIXOUTVEC(config::Info,"Applied field:",applied[0],applied[1],applied[2]);
        FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
        rdt=dt*mat::gamma;
		FIXOUT(config::Info,"Reduced timestep:" << rdt << std::endl);

        mat::sigma = sqrt(2.0*1.38e-23*mat::lambda/(mat::mu*mat::muB*dt*mat::gamma));
        FIXOUT(config::Info,"Sigma prefactor:" << mat::sigma << std::endl);
        llgpf = -1./(1.0+mat::lambda*mat::lambda);
		FIXOUT(config::Info,"Prefactor to LLG equation:" << llgpf << std::endl);

	}
    void integrate(unsigned int& t)
    {
        #ifdef CUDA
        cullg::llgGPU(t);
		if(t%spins::update==0 && rscf)
		{
			spins::calcRealSpaceCorrelationFunction(t);
		}
        #else
        llgCPU::llgCPU(t);
		if(t%spins::update==0 && rscf)
		{
			spins::calcRealSpaceCorrelationFunction(t);
		}
        #endif
    }
}
