// File: fmr.cpp
// Author: Tom Ostler
// Created: 14 June 2013
// Last-modified: 17 Jun 2013 09:10:00
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
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
#include "../inc/exch.h"
void sim::fmr(int argc,char *argv[])
{
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**FMR details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("fmr");
    double redfreq=0.0,Bdrive=0.0;
    int Ncycles=0,MinCycles=0,MaxCycles=0;
    double temp=0.0,ConvVar=2e-11,ConvMean=1e-8;

    std::string opfs;
    try
    {
        setting.lookupValue("temperature",temp);
        setting.lookupValue("redfreq",redfreq);
        setting.lookupValue("Bdrive",Bdrive);
        setting.lookupValue("ConvVar",ConvVar);
        setting.lookupValue("ConvMean",ConvMean);
        setting.lookupValue("Ncycles",Ncycles);
        setting.lookupValue("MaxCycles",MaxCycles);
        setting.lookupValue("MinCycles",MinCycles);
        setting.lookupValue("MagFile",opfs);
    }
    catch(const libconfig::SettingTypeException &stex)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting type error");
    }
    if(redfreq < 1e-12)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Frequency is too small, make sure it is read in properly");
    }

    FIXOUT(config::Info,"Temperature:" << temp << std::endl);
    FIXOUT(config::Info,"Driving frequency/gamma:" << redfreq << std::endl);
    FIXOUT(config::Info,"Magnitude of driving field:" << Bdrive << " [T]" << std::endl);
    FIXOUT(config::Info,"Number of equilibration cycles:" << Ncycles << std::endl);
    FIXOUT(config::Info,"Minimum cycles to average over:" << MinCycles << std::endl);
    FIXOUT(config::Info,"Maximum cycles to average over:" << MaxCycles << std::endl);
    FIXOUT(config::Info,"Average tolerance:" << ConvMean << std::endl);
    FIXOUT(config::Info,"Variance tolerance:" << ConvVar << std::endl);
    FIXOUT(config::Info,"Outputting magnetization data to:" << opfs << std::endl);
    std::ofstream emagfile(opfs.c_str());
    if(!emagfile.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file for writing magnetization data");
    }
    unsigned int ntspc=0;//number of timesteps per cycles
    //----------------------find adjusted frequency----------------------------------------------
    {//only want these variables to have temporary scope
        config::Info << std::setprecision(10);
        double Tp=2.0*M_PI/(redfreq*mat::gamma);//time period
        FIXOUT(config::Info,"Time period of driving field (ideal):" << Tp << std::endl);
        double nts=Tp/llg::dt;//the correct number of timesteps may not be
        FIXOUT(config::Info,"Ideal number of timesteps per cycle:" << nts << std::endl);
        //an integer
        ntspc=int(nts+0.5);//round to nearest integer
        FIXOUT(config::Info,"Actual number of timesteps per cycle:" << ntspc << std::endl);
        double newfreq=2.0*M_PI/(double(ntspc)*llg::dt);
        FIXOUT(config::Info,"New frequency with corrected timesteps:" << newfreq << " rad/s" << std::endl);
        newfreq=newfreq/mat::gamma;
        redfreq=newfreq;
        FIXOUT(config::Info,"New reduced frequency:" << redfreq << std::endl);
    }
    util::RunningStat PC;
    util::TrapezInt TI;
    TI.Clear();
    //TI.bma(2.0*PI*double(cycleaverage)/redfreq);
    TI.bma(double(ntspc)*llg::dt*double(Ncycles));
    PC.Clear();

    FIXOUT(config::Info,"Setting temperature:" << std::flush);
    llg::T=temp;
    SUCCESS(config::Info);
    const double HAppOrig[3]={llg::applied[0],llg::applied[1],llg::applied[2]};
    unsigned int counter=0;
    //Dump the spin configuration when the program exits
    atexit(spins::dump_spins);
    //equilibration
    for(unsigned int cc = 0 ; cc < MinCycles ; cc++)//cycle of driving field loop
    {

        for(unsigned int i = 0 ; i < ntspc ; i++)
        {
            double vf=Bdrive*cos(redfreq*llg::rdt*double(i));
            //std::cout << redfreq << "\t" << llg::rdt << "\t" << double(i) << std::endl;
            llg::applied[0]=HAppOrig[0]+vf;
            llg::applied[1]=HAppOrig[1];
            llg::applied[2]=HAppOrig[2];

            //apply rotation matrix in clockwise direction
            //	[	cos(theta)	sin(theta)	]
            //	[	-sin(theta)	cos(theta)	]
            //  The matrix then takes the vector v and rotates it by an angle theta

            //			v    v'
            //	|     /     *
            //	|    /    *
            //	|   /   *
            //	|  /  *
            //	| / *
            //	|/*_________________________

            //This part was copied from the LLB code for FMR and
            //needs changing if you want to include a rotation of
            //the driving and static magnetic field
            /*double hxt=fields::fH[0];
            double hzt=fields::fH[2];
            fields::fH[0]=hxt*cos(fieldrot)+hzt*sin(fieldrot);
            fields::fH[2]=-hxt*sin(fieldrot)+hzt*cos(fieldrot);
            #endif
            sim::_simt_=counter;
            */
            llg::integrate(counter);
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);

            //calculation of m . dB/dt for FMR integral taking into
            //account rotation in x-z plane
            emagfile << double(counter)*llg::dt/1e-12 << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << "\t" << llg::applied[0] << "\t" << llg::applied[1] << "\t" << llg::applied[2] << std::endl;
            counter++;

        }
    }
    double oldpow=0.0;
    bool broken=false;
    for(unsigned int cc = 0 ; cc < MaxCycles ; cc++)//cycle of driving field loop
    {

        for(unsigned int i = 0 ; i < ntspc ; i++)
        {
            double vf=Bdrive*cos(redfreq*llg::rdt*double(i));
            llg::applied[0]=HAppOrig[0]+vf;
            llg::applied[1]=HAppOrig[1];
            llg::applied[2]=HAppOrig[2];
            llg::integrate(counter);

            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);

            //calculation of m . dB/dt for FMR integral taking into
            //account rotation in x-z plane
            double mdbd=(Bdrive*redfreq*mx*sin(redfreq*llg::rdt*double(i)));
            TI.Push(mdbd);
            emagfile << double(counter)*llg::dt/1e-12 << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << "\t" << llg::applied[0] << "\t" << llg::applied[1] << "\t" << llg::applied[2] << std::endl;
            counter++;
        }


        if(cc%Ncycles==0)
        {
            PC.Push(TI.FinishAndReturn()/(double(ntspc)*llg::dt*double(Ncycles)));
            config::Info << std::setprecision(16);
            config::Info.width(15);config::Info << std::left <<  " | Mean Power = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << PC.Mean();
            config::Info.width(15);config::Info << std::left <<  " | delta P = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(PC.Mean()-oldpow) << " [ " << ConvMean << " ]";
            config::Info.width(15);config::Info << std::left <<  " | Variance = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << PC.Variance() << " [ " << ConvVar << " ] |" << std::endl;
            if((PC.NumDataValues() > 10) && (fabs(PC.Mean()-oldpow) < ConvMean) && (PC.Variance() < ConvVar) && cc > MinCycles)
            {
                broken=true;
                break;
            }
            oldpow=PC.Mean();
            TI.Clear();
        }
    }

}
