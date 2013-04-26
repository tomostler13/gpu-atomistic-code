// File: suscep.h
// Author: Tom Ostler
// Created: 25 Jan 2013
// Last-modified: 26 Apr 2013 12:26:29
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
    setting.lookupValue("MaxRunTime",met);
    setting.lookupValue("EquilTime",et);
    FIXOUT(config::Info,"Maximum run time per temperature:" << met << " seconds" << std::endl);
    unsigned int mrts=int(met/llg::dt),ets=int(et/llg::dt);
    FIXOUT(config::Info,"Maximum run timesteps:" << mrts << std::endl);
    std::string opf;
    setting.lookupValue("SuscepFile",opf);
    FIXOUT(config::Info,"Outputting susceptibility data to:" << opf << std::endl);
    bool outmag=false;
    std::ofstream magout,emagout;
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
    //temperature loop
    for(double T = lT ; T < uT ; T+=dT)
    {
        FIXOUT(config::Info,"Converging temperature:" << T << std::endl);
        llg::T=T;
        if(outmag)
        {
            std::stringstream MagFilesstr,eqsstr;
            MagFilesstr << MagFilestr << "_" << T << ".dat";
            eqsstr << MagFilestr << "_eq_" << T << ".dat";
            std::string mopf=MagFilesstr.str(),eopf=eqsstr.str();
            magout.open(mopf.c_str());
            emagout.open(eopf.c_str());
            if(!magout.is_open() || !emagout.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open file for outputting magnetization data");
            }
        }

        for(unsigned int t = 0 ; t < ets ; t++)
        {
            llg::integrate(t);
            if(t%spins::update==0)
            {
                const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
                const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
                const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
                const double modm=sqrt(mx*mx+my*my+mz*mz);
                emagout << T << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            }
        }
        emagout.close();
        if(emagout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close equilibrium magnetization file");
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

                magout << T << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            }
        }

    }
    std::ofstream fsc("final_spin_config.dat");
    if(!fsc.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not open file for outputting final spin config, outputting to cout");
        for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                {
                    int an=geom::coords(i,j,k,0);
                    if(an>=0)
                    {
                        std::cout << i << "\t" << j << "\t" << k << "\t" << spins::Sx[an] << "\t" << spins::Sy[an] << "\t" << spins::Sz[an] << std::endl;
                    }
                }
            }
        }

    }
    else
    {
        for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                {
                    int an=geom::coords(i,j,k,0);
                    if(an>=0)
                    {
                        fsc << i << "\t" << j << "\t" << k << "\t" << spins::Sx[an] << "\t" << spins::Sy[an] << "\t" << spins::Sz[an] << std::endl;
                    }
                }
            }
        }

        fsc.close();
        if(fsc.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close file for outputting final spin config");
        }
    }


    ofs.close();
    if(ofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close file for outputting magnetization data.");
    }
}

