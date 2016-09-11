// File: llg.cpp
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 11 Sep 2016 16:21:00
#include "../inc/llg.h"
#include "../inc/llgCPU.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/geom.h"
#include "../inc/spins.h"
#include "../inc/util.h"
#include "../inc/rscf.h"
#include <cmath>
#include <sstream>
#include <string>
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace llg
{
    double applied[3]={0,0,0},T,dt,rdt,gyro=1.76e11,muB=9.27e-24,kB=1.38e-23;
    Array<double> Ts,dps,cps;
    Array<double> llgpf;
    std::string scm;

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
        bool llgset=false;
        if(!config::cfg.exists("llg"))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Setting class llg does not exists, setting defaults.");
            llgset=false;
            dt=0.1e-15;
            FIXOUT(config::Info,"Timestep (default):" << dt << " seconds" << std::endl);
            rdt=dt*gyro;
            FIXOUT(config::Info,"Reduced timestep (default):" << rdt << std::endl);
            applied[0]=0.0;applied[1]=0.0;applied[2]=0.0;
            FIXOUTVEC(config::Info,"Applied field (default):",applied[0],applied[1],applied[2]);
            spins::update=100;
            FIXOUT(config::Info,"Update of spins (default):" << spins::update << " [Timesteps]" << std::endl;);
            spins::mag_calc_method=0;
            FIXOUT(config::Info,"Magnetization calculatioin method (default):" << spins::mag_calc_method << std::endl);
            scm="file";
            spins::setSpinsConfig();
            FIXOUT(config::Info,"Initial spin configuratio (default):" << scm << std::endl);
        }
        else
        {
            libconfig::Setting &setting = config::cfg.lookup("llg");
            if(!setting.lookupValue("dt",dt))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Time-step not set (llg.dt (double)), setting default 1e-16 [s]");
                dt=0.1e-15;
                FIXOUT(config::Info,"Timestep (default):" << dt << " seconds" << std::endl);
                rdt=dt*gyro;
                FIXOUT(config::Info,"Reduced timestep (default):" << rdt << std::endl);

            }
            else
            {
                FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
                rdt=dt*gyro;
                FIXOUT(config::Info,"Reduced timestep:" << rdt << std::endl);
            }
            int count=3;
            for(unsigned int i = 0 ; i < 3 ;i++)
            {
                try
                {
                    applied[i]=setting["applied"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errWarnPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " component " << i;
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                    count--;
                }
            }
            if(count==0)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Field not specified, setting default to zero.");
                applied[0]=0.0;applied[1]=0.0;applied[2]=0.0;
                FIXOUTVEC(config::Info,"Applied field (default):",applied[0],applied[1],applied[2]);
            }
            else if(count==3)
            {
                 FIXOUTVEC(config::Info,"Applied field (default):",applied[0],applied[1],applied[2]);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Part of the field was not specified properly, check the warnings.");
            }

            if(!setting.lookupValue("update",spins::update))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read the spin update, setting default to 100 (timesteps)");
                spins::update=100;
                FIXOUT(config::Info,"Update of spins (default):" << spins::update << " [Timesteps]" << std::endl);
            }
            else
            {
                FIXOUT(config::Info,"Update of spins (default):" << spins::update << " [Timesteps]" << std::endl);
            }
            if(!setting.lookupValue("MagnetizationCalculationMethod:",spins::mag_calc_method))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read method for calculating magnetization (llg:MagnetizationCalculationMethod). Setting default to 0 (overall magnetization).");
                spins::mag_calc_method=0;
                FIXOUT(config::Info,"Magnetization calculatioin method (default):" << spins::mag_calc_method << std::endl);

            }
            else
            {
                FIXOUT(config::Info,"Magnetization calculatioin method:" << spins::mag_calc_method << std::endl);
            }
            if(spins::mag_calc_method==9 && geom::Nkset==true)
            {
                util::magDiscSize.resize(3);
                util::magDiscSize.IFill(0);
                for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                {
                    try
                    {
                        util::magDiscSize[xyz]=setting["NoMagDisc"][xyz];
                    }
                    catch(const libconfig::SettingNotFoundException &snf)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream errsstr;
                        errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " component " << xyz;
                        std::string errstr=errsstr.str();
                        error::errMessage(errstr);
                    }
                }
                FIXOUTVEC(config::Info,"Number of intervals for calculating discrete magnetization:",util::magDiscSize[0],util::magDiscSize[1],util::magDiscSize[2]);

                double fractpart,intpart;

                double no_kpx_pd=static_cast<double>(geom::dim[0]*geom::Nk[0])/static_cast<double>(util::magDiscSize[0]);
                fractpart=std::modf(no_kpx_pd,&intpart);
                if(fabs(fractpart)>1e-12)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    printf ("No. points per discretization cell %f = %f + %f (the fractional part is bigger than 1e-12) \n", no_kpx_pd, intpart, fractpart);
                    error::errMessage("Error in discretization (in x). At the moment the code requires a perfect number of points so that geom::dim*geom::Nk/magDiscSize is EXACTLY an integer (to an accuracy of 1e-12).");
                }
                double no_kpy_pd=static_cast<double>(geom::dim[1]*geom::Nk[1])/static_cast<double>(util::magDiscSize[1]);
                fractpart=std::modf(no_kpy_pd,&intpart);
                if(fabs(fractpart)>1e-12)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    printf ("No. points per discretization cell %f = %f + %f (the fractional part is bigger than 1e-12) \n", no_kpy_pd, intpart, fractpart);
                    error::errMessage("Error in discretization (in y). At the moment the code requires a perfect number of points so that geom::dim*geom::Nk/magDiscSize is EXACTLY an integer (to an accuracy of 1e-12).");
                }
                double no_kpz_pd=static_cast<double>(geom::dim[2]*geom::Nk[2])/static_cast<double>(util::magDiscSize[2]);
                fractpart=std::modf(no_kpz_pd,&intpart);
                if(fabs(fractpart)>1e-12)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Error in discretization (in z). At the moment the code requires a perfect number of points so that geom::dim*geom::Nk/magDiscSize is EXACTLY an integer (to an accuracy of 1e-12).");
                    printf ("No. points per discretization cell %f = %f + %f (the fractional part is bigger than 1e-12) \n", no_kpz_pd, intpart, fractpart);
                }
                if(setting.lookupValue("DiscOutputFormat",util::disOutForm))
                {
                    FIXOUT(config::Info,"Output format for mag disc output (MagnetizationCalculationMethod):" << util::disOutForm << std::endl);
                }
                else
                {
                    error::errWarnPreamble(__FILE__,__LINE__);
                    error::errWarning("Could not read the DiscOutputFormat (llg:DiscOutputFormat (string)), defaulting to VTU.");
                    util::disOutForm="vtu";
                }
            }
            else if(spins::mag_calc_method==9 && geom::Nkset==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("To use magnetisation method 9 (discretization into cells) you have to define Nm.");
            }
            if(!setting.lookupValue("OutputMagnetization",spins::output_mag))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Output overall magnetization (llg.OutputMagnetization) not specified. Defaulting to true.");
            }
            FIXOUT(config::Info,"Initializing output of magnetization:" << std::endl);
            if(setting.lookupValue("SpinConfigMethod",scm))
            {
                FIXOUT(config::Info,"How are spins initially configured?:" << scm << std::endl);
                if(scm=="random")
                {
                    spins::setSpinsRandom();
                }
                else if(scm=="file")
                {
                    spins::setSpinsConfig();
                }
                else if(scm=="chequerx")
                {
                    spins::setSpinsChequerX();
                }
                else if(scm=="chequery")
                {
                    spins::setSpinsChequerY();
                }
                else if(scm=="chequerz")
                {
                    spins::setSpinsChequerZ();
                }
                else if(scm=="species")
                {
                    spins::setSpinsSpecies(setting);
                }
                else if(scm=="vampire")
                {
                    spins::setSpinsVampire(setting);
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Method for initialising spins not recognised (llg.SpinConfigMethod (string). Check you config file.");
                }
            }
            else
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read the method for applying the initial spin config. Setting llg.SpinConfigMethod (string), defaulting to file.");
                scm="file";
                spins::setSpinsConfig();
                FIXOUT(config::Info,"Initial spin configuratio (default):" << scm << std::endl);
            }
        }
        util::init_output();

        if(geom::ucm.NumAtomsUnitCell() > 5)
        {
            config::openLogFile();
        }

        config::printline(config::Info);
        //Output the prefactor to the LLG for each species to the output file
        for(unsigned int i = 0 ; i < geom::ucm.NumAtomsUnitCell() ; i++)
        {
            std::stringstream sstr,sstr1;
            std::string str=sstr.str();
            sstr << "Sigma prefactor for unit cell atom " << i << ":";
            sstr1 << "Prefactor for LLG for unit cell atom " << i << ":";
            geom::ucm.SetSigma(i,sqrt(2.0*kB*geom::ucm.GetDamping(i)/(geom::ucm.GetMu(i)*muB*dt*geom::ucm.GetGamma(i)*gyro)));
            geom::ucm.Setllgpf(i,-1./(1.0+geom::ucm.GetDamping(i)*geom::ucm.GetDamping(i)));

            if(i < 5)
            {
                FIXOUT(config::Info,str.c_str() << geom::ucm.GetSigma(i) << std::endl);
                str=sstr1.str();
                FIXOUT(config::Info,str.c_str() << geom::ucm.Getllgpf(i) << std::endl);
                config::printline(config::Info);
            }
            if(geom::ucm.NumAtomsUnitCell() > 5 && geom::logunit)
            {
                FIXOUT(config::Log,str.c_str() << geom::ucm.GetSigma(i) << std::endl);
                str=sstr1.str();
                FIXOUT(config::Log,str.c_str() << geom::ucm.Getllgpf(i) << std::endl);
                config::printline(config::Log);
            }
        }
        if(geom::ucm.NumAtomsUnitCell()>5)
        {
            for(unsigned int i = 0 ; i < 4 ; i++)
            {
                FIXOUT(config::Info,"             . . . " << "    . . ." << std::endl);
            }
            if(geom::logunit)
            {
                FIXOUT(config::Info,"FOR COMPLETE LLG INFO FOR UNIT CELL SEE LOG FILE:" << "   log.dat" << std::endl);
            }
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
            geom::sigma[i]=geom::ucm.GetSigma(aiuc);
        }
        llg::Ts.resize(geom::ucm.GetNMS());
        llg::dps.resize(geom::ucm.GetNMS());
        llg::cps.resize(geom::ucm.GetNMS());

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
