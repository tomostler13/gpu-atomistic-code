// File: llg.cpp
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 25 Sep 2014 10:42:38
#include "../inc/llg.h"
#include "../inc/llgCPU.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/geom.h"
#include <cmath>
#include <sstream>
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace llg
{
    double applied[3]={0,0,0},T,dt,rdt,gyro=1.76e11,muB=9.27e-24;
    Array<double> llgpf;

	void initLLG(int argc,char *argv[])
	{
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**LLG details***" << std::endl;

        //resize the llgpf array
        llgpf.resize(geom::nms);
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
        for(unsigned int i = 0 ; i < 3 ;i++)
        {
            applied[i]=setting["applied"][i];
        }
        FIXOUTVEC(config::Info,"Applied field:",applied[0],applied[1],applied[2]);
        FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
        rdt=dt*gyro;
		FIXOUT(config::Info,"Reduced timestep:" << rdt << std::endl);
        //set the prefactor of the LLG for each species
        for(unsigned int i = 0 ; i < geom::nms ; i++)
        {
            //mat::sigma[i] = sqrt(2.0*1.38e-23*mat::lambda[i]/(mat::mu[i]*mat::muB*dt*mat::gyro*mat::gamma[i]));
            std::stringstream sstr;
            sstr << "Sigma prefactor for species " << i << ":";
            std::string str=sstr.str();
            //FIXOUT(config::Info,str.c_str() << mat::sigma[i] << std::endl);
            //llgpf[i] = -1./(1.0+mat::lambda[i]*mat::lambda[i]);
            sstr.str("");
            sstr << "Prefactor for LLG for species " << i << ":";
            str=sstr.str();
            FIXOUT(config::Info,str.c_str() << llgpf[i] << std::endl);
        }

	}
    void integrate(unsigned int& t)
    {
        #ifdef CUDA
        cullg::llgGPU(t);
        #else
        llgCPU::llgCPU(t);
        #endif
    }
}
