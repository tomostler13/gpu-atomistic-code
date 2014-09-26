// File: llg.cpp
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 26 Sep 2014 11:39:08
#include "../inc/llg.h"
#include "../inc/llgCPU.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/geom.h"
#include "../inc/spins.h"
#include <cmath>
#include <sstream>
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace llg
{
    double applied[3]={0,0,0},T,dt,rdt,gyro=1.76e11,muB=9.27e-24,kB=1.38e-23;
    Array<double> llgpf;

	void initLLG(int argc,char *argv[])
	{
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**LLG details***" << std::endl;

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
        setting.lookupValue("update",spins::update);
        FIXOUT(config::Info,"Spin update:" << spins::update << " [Timesteps]" << std::endl);
        FIXOUTVEC(config::Info,"Applied field:",applied[0],applied[1],applied[2]);
        FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
        rdt=dt*gyro;
		FIXOUT(config::Info,"Reduced timestep:" << rdt << std::endl);
        setting.lookupValue("MagnetizationCalculationMethod:",spins::mag_calc_method);
        if(geom::ucm.NumAtomsUnitCell() > 5)
        {
            config::openLogFile();
        }

        config::printline(config::Info);
        //Output the prefactor to the LLG for each species to the output file
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            std::stringstream sstr,sstr1;
            sstr << "Sigma prefactor for unit cell atom " << i << ":";
            std::string str=sstr.str();
            geom::ucm.SetSigma(i,sqrt(2.0*kB*geom::ucm.GetDamping(i)/(geom::ucm.GetMu(i)*muB*dt*geom::ucm.GetGamma(i)*gyro)));
            sstr1 << "Prefactor for LLG for unit cell atom " << i << ":";
            geom::ucm.Setllgpf(i,-1./(1.0+geom::ucm.GetDamping(i)*geom::ucm.GetDamping(i)));

            if(i < 5)
            {
                FIXOUT(config::Info,str.c_str() << geom::ucm.GetSigma(i) << std::endl);
                FIXOUT(config::Info,str.c_str() << geom::ucm.Getllgpf(i) << std::endl);
                config::printline(config::Info);
                str=sstr1.str();
            }
            if(geom::ucm.NumAtomsUnitCell() > 5)
            {
                FIXOUT(config::Log,str.c_str() << geom::ucm.GetSigma(i) << std::endl);
                FIXOUT(config::Info,str.c_str() << geom::ucm.Getllgpf(i) << std::endl);
                config::printline(config::Log);
            }
        }
        if(geom::ucm.NumAtomsUnitCell()>5)
        {
            for(unsigned int i = 0 ; i < 4 ; i++)
            {
                FIXOUT(config::Info,"             . . . " << "    . . ." << std::endl);
            }
            FIXOUT(config::Info,"FOR COMPLETE LLG INFO FOR UNIT CELL SEE LOG FILE:" << "   log.dat" << std::endl);
            for(unsigned int i = 0 ; i < 4 ; i++)
            {
                FIXOUT(config::Info,"             . . . " << "    . . ." << std::endl);
            }

        }
        FIXOUT(config::Info,"Setting the llg and thermal prefactors to 1D arrays:" << std::flush);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int aiuc=geom::lu(i,4);;
            geom::llgpf[i]=geom::ucm.Getllgpf(aiuc);
            geom::sigma[i]=geom::ucm.Getllgpf(aiuc);
        }
        SUCCESS(config::Info);

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
