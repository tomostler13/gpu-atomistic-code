// File: mvt.h
// Author: Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 22 Jun 2016 14:06:42
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "../inc/arrays.h"
#include "../inc/defines.h"
#include "../inc/util.h"
#include "../inc/config.h"
#include "../inc/spins.h"
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
    util::RunningStat *MS=new util::RunningStat[geom::ucm.GetNMS()];
    util::ofs << "#The magnetization is output into indices for each temperature and depending on the output method" << std::endl;
    if(dT>0.0)
    {
        //temperature loop
        for(double T = lT ; T < uT ; T+=dT)
        {
            double modm[geom::ucm.GetNMS()],oldmean[geom::ucm.GetNMS()];
            bool convTF[geom::ucm.GetNMS()];
            config::printline(config::Info);
            FIXOUT(config::Info,"Converging temperature:" << T << std::endl);
            llg::T=T;
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                MS[i].Clear();
                modm[i]=0.0;
                oldmean[i]=0.0;
                convTF[i]=false;
            }
            for(unsigned int t = 0 ; t < eminrts ; t++)
            {
                llg::integrate(t);
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);
                    util::outputSpinsVTU(t);
                }
            }
            //have all the magnetization converged?
            bool allconv=false;
            for(unsigned int t = eminrts ; t < eminrts+mrts ; t++)
            {
                llg::integrate(t);
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);

                    /*                if(t>int(10e-12/llg::dt))
                                      {*/

                    unsigned int counter=0;
                    for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
                    {
                        modm[i]=sqrt(spins::mag(i,0)*spins::mag(i,0)+spins::mag(i,1)*spins::mag(i,1)+spins::mag(i,2)*spins::mag(i,2));
                        MS[i].Push(modm[i]);
                        config::Info.width(15);config::Info << "| Species " << i << " mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS[i].Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MS[i].Mean()-oldmean[i]) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS[i].Variance() << " [ " << convvar << " ]|" << std::endl;
                        if(((fabs(MS[i].Mean()-oldmean[i])) < convmean) && (MS[i].Variance()<convvar))
                        {
                            convTF[i]=true;
                            counter++;
                        }
                        oldmean[i]=MS[i].Mean();
                    }
                    if(counter==geom::ucm.GetNMS())
                    {
                        allconv=true;
                        break;
                    }

                    //}

                }
            }
            ofs << T;
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                ofs << "\t" << modm[i] << std::flush;
            }
            ofs << std::endl;

            FIXOUT(config::Info,"Converged?" << config::isTF(allconv) << std::endl);
            util::ofs << std::endl << std::endl;
        }
    }
    if(dT<0.0)
    {
        //temperature loop
        for(double T = uT ; T > lT ; T+=dT)
        {
            double modm[geom::ucm.GetNMS()],oldmean[geom::ucm.GetNMS()];
            bool convTF[geom::ucm.GetNMS()];
            config::printline(config::Info);
            FIXOUT(config::Info,"Converging temperature:" << T << std::endl);
            llg::T=T;
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                MS[i].Clear();
                modm[i]=0.0;
                oldmean[i]=0.0;
                convTF[i]=false;
            }
            for(unsigned int t = 0 ; t < eminrts ; t++)
            {
                llg::integrate(t);
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);
                }
            }
            //have all the magnetization converged?
            bool allconv=false;
            for(unsigned int t = eminrts ; t < eminrts+mrts ; t++)
            {
                llg::integrate(t);
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);

                    /*                if(t>int(10e-12/llg::dt))
                                      {*/

                    unsigned int counter=0;
                    for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
                    {
                        modm[i]=sqrt(spins::mag(i,0)*spins::mag(i,0)+spins::mag(i,1)*spins::mag(i,1)+spins::mag(i,2)*spins::mag(i,2));
                        MS[i].Push(modm[i]);
                        config::Info.width(15);config::Info << "| Species " << i << " mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS[i].Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MS[i].Mean()-oldmean[i]) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MS[i].Variance() << " [ " << convvar << " ]|" << std::endl;
                        if(((fabs(MS[i].Mean()-oldmean[i])) < convmean) && (MS[i].Variance()<convvar))
                        {
                            convTF[i]=true;
                            counter++;
                        }
                        oldmean[i]=MS[i].Mean();
                    }
                    if(counter==geom::ucm.GetNMS())
                    {
                        allconv=true;
                        break;
                    }

                    //}

                }
            }
            ofs << T;
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                ofs << "\t" << modm[i] << std::flush;
            }
            ofs << std::endl;

            FIXOUT(config::Info,"Converged?" << config::isTF(allconv) << std::endl);
            util::ofs << std::endl << std::endl;
        }
    }
    ofs.close();
    if(ofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close file for outputting magnetization data.");
    }
}
