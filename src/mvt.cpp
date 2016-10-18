// File: mvt.h
// Author: Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 18 Oct 2016 13:16:40
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
void sim::MvT()
{
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**M(T) details***" << std::endl;

    if(!config::cfg.exists("mvt"))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting \"mvt\" does not exist. Check your config file.");
    }
    libconfig::Setting &setting = config::cfg.lookup("mvt");
    double lT=0.0,uT=0.0,dT=0.0,convmean=0.0,convvar=0.0,met=0.0,minrt=0.0;
    if(!setting.lookupValue("Lower_temp",lT))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the lower temperature mvt.Lower_temp (double). It is required.");
    }
    else
    {
        FIXOUT(config::Info,"Lower temperature:" << lT << std::endl);
    }
    if(!setting.lookupValue("Upper_temp",uT))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the upper temperature (mvt.Upper_temp (double)). It is required");
    }
    else
    {
        FIXOUT(config::Info,"Upper temperature:" << uT << std::endl);
    }
    if(!setting.lookupValue("Temp_step",dT))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the step in temperature (mvt.Temp_step (double)). It is required");
    }
    else
    {
        FIXOUT(config::Info,"Temperature step:" << dT << std::endl);
    }
    if(!setting.lookupValue("Mean_tol",convmean))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Mean tolerance could not be read (mvt.Mean_tol (double)), defaulting to 1e-5");
        convmean=1e-5;
        FIXOUT(config::Info,"Converging mean to (default):" << convmean << std::endl);
    }
    else
    {
        FIXOUT(config::Info,"Converging mean to:" << convmean << std::endl);
    }
    if(!setting.lookupValue("Variance_tol",convvar))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Variance tolerance could not be read (mvt.Variance_tol (double)), defaulting to 1e-6");
        convvar=1e-6;
    }
    else
    {
        FIXOUT(config::Info,"Converging variance to:" << convvar << std::endl);
    }
    if(!setting.lookupValue("MaxRunTime",met))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Maximum run time per temperature step not read (mvt.MaxRunTime (double)), defaulting to 25ps");
        met=25e-12;
        FIXOUT(config::Info,"Maximum run time per temperature (default):" << met << " seconds" << std::endl);
    }
    else
    {
        FIXOUT(config::Info,"Maximum run time per temperature:" << met << " seconds" << std::endl);
    }
    if(!setting.lookupValue("MinRunTime",minrt))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Minimum run time per temperature step not read (mvt.MinRunTime (double)), defaulting to 5ps");
        met=5e-12;
        FIXOUT(config::Info,"Minimum run time per temperature (default):" << met << " seconds" << std::endl);
    }
    else
    {
        FIXOUT(config::Info,"Minimum runtime per temperature step:" << minrt << std::endl);
    }
    unsigned int mrts=int(met/llg::dt),eminrts=int(minrt/llg::dt);
    FIXOUT(config::Info,"Maximum run timesteps:" << mrts << std::endl);
    FIXOUT(config::Info,"Min rum timesteps:" << eminrts << std::endl);
    std::string opf;
    if(!setting.lookupValue("MvTFile",opf))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not read the name of the outputfile for the magnetization. Defaulting to MvT.dat");
        opf="MvT.dat";
        FIXOUT(config::Info,"Outputting magnetization data to (default):" << opf << std::endl);
    }
    else
    {
        FIXOUT(config::Info,"Outputting magnetization data to (default):" << opf << std::endl);
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

    //Output VTU files? Two options exist 1) regularly outputting the spin configuration at given intervals, or 2) outputting at the end for each temperature. The default is 2.
    bool outVTU=false;
    unsigned int outVTUtype=2,outVTUfreq=1;
    if(!setting.lookupValue("OutputSpinsVTU",outVTU))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not read if you want to output the VTU files during M(T) calculation (mvt.OutputSpinsVTU (bool)), defaulting to false.");
        outVTU=false;
        FIXOUT(config::Info,"Outputting spins in VTU format (default)?:" << config::isTF(outVTU) << std::endl);
    }
    else
    {
        FIXOUT(config::Info,"Outputting spins in VTU format?:" << config::isTF(outVTU) << std::endl);
        if(!setting.lookupValue("OutputVTUType",outVTUtype))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Could not read output type of VTUs, (mvt.OutputVTUType (int)), defaulting to 2 (at the end of each T step)");
            outVTUtype=2;
            FIXOUT(config::Info,"Output VTU type (default):" << outVTUtype << std::endl);
        }
        else
        {
            FIXOUT(config::Info,"Output VTU type:" << outVTUtype << std::endl);
        }
        if(outVTUtype==1)
        {
            if(!setting.lookupValue("OutputVTUFreq",outVTUfreq))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read output VTU frequency (mvt.OutputVTUFreq (int) in units of update, defaulting to 100");
                outVTUfreq=100;
                FIXOUT(config::Info,"Output frequency of VTUs (default):" << outVTUfreq << std::endl);
            }
        }
        else
        {
            FIXOUT(config::Info,"Output frequency of VTUs:" << outVTUfreq << std::endl);
        }
    }

    util::RunningStat *MS=new util::RunningStat[geom::ucm.GetNMS()];
    util::ofs << "#The magnetization is output into indices for each temperature and depending on the output method" << std::endl;
    unsigned int globt=0.0;
    unsigned int VTUcount=0;
    if(dT>0.0)
    {
        //temperature loop
        for(double T = lT ; T < uT ; T+=dT)
        {
            double *modm=NULL,*oldmean=NULL;
            modm = new double[geom::ucm.GetNMS()];
            oldmean = new double[geom::ucm.GetNMS()];
            bool *convTF=NULL;
            convTF = new bool[geom::ucm.GetNMS()];
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
                globt++;
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);
                    if(outVTUtype==1 && VTUcount==outVTUfreq)
                    {
                        util::outputSpinsVTU(t,T);
                        VTUcount=0;
                    }
                    else if(outVTUtype==1 && VTUcount!=outVTUfreq)
                    {
                        VTUcount++;
                    }
                }
            }
            //have all the magnetization converged?
            bool allconv=false;
            for(unsigned int t = eminrts ; t < eminrts+mrts ; t++)
            {
                llg::integrate(t);
                globt++;
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);
                    if(outVTUtype==1 && VTUcount==outVTUfreq)
                    {
                        util::outputSpinsVTU(t,T);
                        VTUcount=0;
                    }
                    else if(outVTUtype==1 && VTUcount!=outVTUfreq)
                    {
                        VTUcount++;
                    }

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
            if(outVTUtype==2)
            {
                util::outputSpinsVTU(T);
            }
            ofs << T;
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                ofs << "\t" << modm[i] << std::flush;
            }
            ofs << std::endl;

            FIXOUT(config::Info,"Converged?" << config::isTF(allconv) << std::endl);
            util::ofs << std::endl << std::endl;
            try
            {
                delete [] oldmean;
                oldmean=NULL;
            }
            catch(...)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not delete array oldmean");
            }
            try
            {
                delete [] convTF;
                convTF=NULL;
            }
            catch(...)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not delete array convTF");
            }
            try
            {
                delete [] modm;
                modm=NULL;
            }
            catch(...)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not delete array modm");
            }
        }
    }
    if(dT<0.0)
    {
        //temperature loop
        for(double T = uT ; T > lT ; T+=dT)
        {
            double *modm=NULL,*oldmean=NULL;
            modm = new double[geom::ucm.GetNMS()];
            oldmean = new double[geom::ucm.GetNMS()];
            bool *convTF=NULL;
            convTF = new bool[geom::ucm.GetNMS()];
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
                globt++;
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);
                    if(outVTUtype==1 && VTUcount==outVTUfreq)
                    {
                        util::outputSpinsVTU(t,T);
                        VTUcount=0;
                    }
                    else if(outVTUtype==1 && VTUcount!=outVTUfreq)
                    {
                        VTUcount++;
                    }
                }
            }
            //have all the magnetization converged?
            bool allconv=false;
            for(unsigned int t = eminrts ; t < eminrts+mrts ; t++)
            {
                llg::integrate(t);
                globt++;
                if(t%spins::update==0)
                {
                    util::calc_mag();
                    util::output_mag(t);
                    util::outputSpinsVTU(globt);
                    if(outVTUtype==1 && VTUcount==outVTUfreq)
                    {
                        util::outputSpinsVTU(t,T);
                        VTUcount=0;
                    }
                    else if(outVTUtype==1 && VTUcount!=outVTUfreq)
                    {
                        VTUcount++;
                    }

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
            if(outVTUtype==2)
            {
                util::outputSpinsVTU(T);
            }
            ofs << T;
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                ofs << "\t" << modm[i] << std::flush;
            }
            ofs << std::endl;

            FIXOUT(config::Info,"Converged?" << config::isTF(allconv) << std::endl);
            util::ofs << std::endl << std::endl;
            try
            {
                delete [] oldmean;
                oldmean=NULL;
            }
            catch(...)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not delete array oldmean");
            }
            try
            {
                delete [] convTF;
                convTF=NULL;
            }
            catch(...)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not delete array convTF");
            }
            try
            {
                delete [] modm;
                modm=NULL;
            }
            catch(...)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not delete array modm");
            }
        }
    }
    ofs.close();
    if(ofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close file for outputting magnetization data.");
    }
}
