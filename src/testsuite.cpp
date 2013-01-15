// File: geom.cpp
// Author:Tom Ostler
// Last-modified: 03 Jan 2013 13:59:01
//------header files for globals-------
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/mat.h"
#include "../inc/intmat.h"
#include "../inc/spins.h"
#include "../inc/fields.h"
#include "../inc/maths.h"
#include "../inc/util.h"
#include "../inc/neigh.h"
#include "../inc/LLBCPU.h"
#include "../inc/LLB.h"
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
#include "../inc/tdp.h"
#include "../inc/sim.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
void sim::testsuite(int argc,char *argv[])
{
    assert(geom::gi);
    assert(tdp::tdpi);
    assert(fields::fi);
    assert(spins::si);
    assert(util::ui);
    assert(llb::LLBi);
    assert(mat::mi);
    assert(neigh::ni);
    assert(config::lcf);
    assert(sim::simi);

    //min temperature
    double MvsTminTemp = 0.0000001;
    //maximum temperature
    double MvsTmaxTemp = 700.0;
    //temperature step
    double MvsTdTemp = 5.0;
    //convergence mean
    double convmean = 0.0;
    //variance to converge to
    double convvar = 0.0;
    //string containing output filename
    std::string MvsTopf;
    std::string MvsTinfo;
    std::string oneDboltzfile;
    std::ofstream MvsTopfs;
    std::ofstream MvsTinfofs;
    std::ofstream bfofs;
    bool MvsTCase1=false;
    bool MvsTCase2=false;
    bool MvsTCase3=false;
    bool MvsTCase4=false;
    bool boltz1d=true;
    unsigned int b1dts=0;
    double b1dtemp=1e-7;


    if(tdp::toh!="uniform")
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Type of heating (toh) is not set to uniform, required for the test suite");
    }

    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Simulation details***" << std::endl;
    FIXOUT(config::Info,"Simulation type:" << sim::sim_type << std::endl);

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

    libconfig::Setting &setting = config::cfg.lookup("sim.testsuite");
    try
    {
        setting.lookupValue("MvsTminTemp",MvsTminTemp);
        setting.lookupValue("MvsTmaxTemp",MvsTmaxTemp);
        setting.lookupValue("MvsTdTemp",MvsTdTemp);
        setting.lookupValue("MvsToutputFile",MvsTopf);
        setting.lookupValue("MvsTInfoFile",MvsTinfo);
        setting.lookupValue("MvsTConvMean",convmean);
        setting.lookupValue("MvsTConvVar",convvar);
        setting.lookupValue("oneDBoltzFile",oneDboltzfile);
        setting.lookupValue("MvsTCase1",MvsTCase1);
        setting.lookupValue("MvsTCase2",MvsTCase2);
        setting.lookupValue("MvsTCase3",MvsTCase3);
        setting.lookupValue("MvsTCase4",MvsTCase4);
        setting.lookupValue("boltz1d",boltz1d);
        setting.lookupValue("boltz1dtimesteps",b1dts);
        setting.lookupValue("boltz1dTemp",b1dtemp);
    }
    catch(const libconfig::SettingTypeException &stex)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting type error");
    }
    FIXOUT(config::Info,"MvsT minimum temperature:" << MvsTminTemp << std::endl);
    FIXOUT(config::Info,"MvsT maximum temperature:" << MvsTmaxTemp << std::endl);
    FIXOUT(config::Info,"MvsT change in temperature:" << MvsTdTemp << std::endl);
    FIXOUT(config::Info,"MvsT output file:" << MvsTopf << std::endl);
    FIXOUT(config::Info,"MvsT info file:" << MvsTinfo << std::endl);
    FIXOUT(config::Info,"MvsT Case 1:" << config::isTF(MvsTCase1) << std::endl);
    FIXOUT(config::Info,"MvsT Case 2:" << config::isTF(MvsTCase2) << std::endl);
    FIXOUT(config::Info,"MvsT Case 3:" << config::isTF(MvsTCase3) << std::endl);
    FIXOUT(config::Info,"MvsT Case 4:" << config::isTF(MvsTCase4) << std::endl);
    FIXOUT(config::Info,"1D Boltzmann test file:" << oneDboltzfile << std::endl);
    FIXOUT(config::Info,"1D Boltzmann test number of steps:" << b1dts << std::endl);
    FIXOUT(config::Info,"1D Boltzmann test temperature:" << b1dtemp << std::endl);
    if(convmean<1e-20)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("MvsT mean convergence far too small");
    }
    if(convvar<1e-20)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("MvsT variance convergence far too small");
    }

    FIXOUT(config::Info,"MvsT mean converged to:" << convmean << std::endl);
    FIXOUT(config::Info,"MvsT variance converged to:" << convvar << std::endl);

    MvsTopfs.open(MvsTopf.c_str());
    MvsTinfofs.open(MvsTinfo.c_str());
    bfofs.open(oneDboltzfile.c_str());
    if(!MvsTopfs.is_open() || !MvsTinfofs.is_open() || !bfofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Cannot open file for writing");
    }
    MvsTopfs << "#The MvsT tests include switch the flags:" << std::endl;
    MvsTopfs << "#-fields::checkdipolar\n#-fields::checkthermal\n#-fields::checkexchange" << std::endl;
    MvsTopfs << "#The MvsT data is output in indexes that can be indexed by name as follows" << std::endl;
    MvsTopfs << "#1 - no dipolar, no thermal, no exchange" << std::endl;
    MvsTopfs << "#2 - no dipolar, no thermal, with exchange" << std::endl;
    MvsTopfs << "#3 - no dipolar, with thermal, no exchange" << std::endl;
    MvsTopfs << "#4 - no dipolar, with thermal, with exchange" << std::endl;
    MvsTopfs << "#T [K]\t<M>\n";
    MvsTopfs << "\n\n" << std::endl;

    //initialize spins
    FIXOUT(MvsTinfofs,"Initializing spins to:" << "(0.1,0.0," << spins::sz[0] << ")" << std::endl);
    for(unsigned int i = 0 ; i < geom::ss ; i++)
    {
        spins::sx[i]=0.1;
        spins::sy[i]=0.0;
        spins::sz[i]=sqrt(1.0-spins::sx[i]*spins::sx[i]);
    }
    //set fields to (0,0,0)
    FIXOUT(MvsTinfofs,"Setting applied field to:" << "(0,0,0)" << std::endl);
    for(unsigned int i = 0 ; i < 3 ; i++)
    {
        fields::fH[i]=0.0;
    }

    // Case 1 - No dipolar, no thermal, no exchange
    config::printline(MvsTinfofs);
    if(MvsTCase1==true)
    {
        MvsTinfofs << "Running MvsT case 1" << std::endl;
        MvsTopfs << "#one" << std::endl;

        fields::checkdipolar=false;
        fields::checkthermal=false;
        fields::checkexchange=false;
        tdp::exchstifffp=tdp::exchstiffNULL;
        MvsTopfs << std::setprecision(16);

        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            spins::sx[i]=0.1;
            spins::sy[i]=0.0;
            spins::sz[i]=sqrt(1.0-spins::sx[i]*spins::sx[i]);
        }

        for(double T = MvsTminTemp ; T < MvsTmaxTemp ; T+=MvsTdTemp)
        {
            //set the uniform system temperature
            tdp::uniformtemp=T;
            util::RunningStat MC;
            MC.Clear();
            MvsTinfofs << std::endl;
            FIXOUT(MvsTinfofs,"Converging temperatue:" << T << std::endl << std::endl);
            double oldmean=0.0;
            #ifdef CUDA
            cullb::initGPU();
            #endif
            for(unsigned int t = 0 ; ; )
            {
                util::calcm();
                MC.Push(spins::modm);
                if(t>0 && t%1000==0)
                {
                    MvsTinfofs.width(15);MvsTinfofs << "| Mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MC.Mean()-oldmean) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Variance() << " [ " << convvar << " ] " << std::endl;
                    if(((fabs(MC.Mean()-oldmean)) < convmean) && (MC.Variance()<convvar))
                    {
                        MvsTopfs << T << "\t" << spins::modm << std::endl;
                        break;
                    }
                    oldmean=MC.Mean();
                }
                LLB::integrate(t);
    //            std::cout << t << "\t" << spins::mx << "\t" << spins::my << "\t" << spins::mz << "\t" << spins::modm << std::endl;

            }
            MvsTinfofs << std::endl;

        }
    }
    if(MvsTCase2==true)
    {
        // Case 2 - No dipolar, no thermal, with exchange
        config::printline(MvsTinfofs);
        MvsTinfofs << "Running MvsT case 2" << std::endl;
        fields::checkdipolar=false;
        fields::checkthermal=false;
        fields::checkexchange=true;
        tdp::exchstifffp=tdp::exchstiff0;
        MvsTopfs << std::setprecision(16);
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            spins::sx[i]=0.1;
            spins::sy[i]=0.0;
            spins::sz[i]=sqrt(1.0-spins::sx[i]*spins::sx[i]);
        }
        MvsTopfs << std::endl << std::endl;
        MvsTopfs << "#two" << std::endl;
        for(double T = MvsTminTemp ; T < MvsTmaxTemp ; T+=MvsTdTemp)
        {
            //set the uniform system temperature
            tdp::uniformtemp=T;
            util::RunningStat MC;
            MC.Clear();
            MvsTinfofs << std::endl;
            FIXOUT(MvsTinfofs,"Converging temperatue:" << T << std::endl << std::endl);
            double oldmean=0.0;
            #ifdef CUDA
            cullb::initGPU();
            #endif

            for(unsigned int t = 0 ; ; )
            {
                util::calcm();
                MC.Push(spins::modm);
                if(t>0 && t%1000==0)
                {
                    MvsTinfofs.width(15);MvsTinfofs << "| Mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MC.Mean()-oldmean) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Variance() << " [ " << convvar << " ] " << std::endl;
                    if(((fabs(MC.Mean()-oldmean)) < convmean) && (MC.Variance()<convvar))
                    {
                        MvsTopfs << T << "\t" << spins::modm << std::endl;
                        break;
                    }
                    oldmean=MC.Mean();
                }
                LLB::integrate(t);

    //            std::cout << t << "\t" << spins::mx << "\t" << spins::my << "\t" << spins::mz << "\t" << spins::modm << std::endl;

            }
            MvsTinfofs << std::endl;

        }

        config::printline(MvsTinfofs);
    }
    if(MvsTCase3==true)
    {
        // Case 3 - No dipolar, with thermal, no exchange
        config::printline(MvsTinfofs);
        MvsTinfofs << "Running MvsT case 3" << std::endl;
        fields::checkdipolar=false;
        fields::checkthermal=true;
        fields::checkexchange=false;
        tdp::exchstifffp=tdp::exchstiffNULL;
        MvsTopfs << std::setprecision(16);
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            spins::sx[i]=0.1;
            spins::sy[i]=0.0;
            spins::sz[i]=sqrt(1.0-spins::sx[i]*spins::sx[i]);
        }
        MvsTopfs << std::endl << std::endl;
        MvsTopfs << "#three" << std::endl;
        for(double T = MvsTminTemp ; T < MvsTmaxTemp ; T+=MvsTdTemp)
        {
            //set the uniform system temperature
            tdp::uniformtemp=T;
            util::RunningStat MC;
            MC.Clear();
            MvsTinfofs << std::endl;
            FIXOUT(MvsTinfofs,"Converging temperatue:" << T << std::endl << std::endl);
            double oldmean=0.0;
            #ifdef CUDA
            cullb::initGPU();
            #endif

            for(unsigned int t = 0 ; ; )
            {
                util::calcm();
                MC.Push(spins::modm);
                if(t>0 && t%1000==0)
                {
                    MvsTinfofs.width(15);MvsTinfofs << "| Mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MC.Mean()-oldmean) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Variance() << " [ " << convvar << " ] " << std::endl;
                    if(((fabs(MC.Mean()-oldmean)) < convmean) && (MC.Variance()<convvar))
                    {
                        MvsTopfs << T << "\t" << spins::modm << std::endl;
                        break;
                    }
                    oldmean=MC.Mean();
                }
                LLB::integrate(t);

    //            std::cout << t << "\t" << spins::mx << "\t" << spins::my << "\t" << spins::mz << "\t" << spins::modm << std::endl;

            }
            MvsTinfofs << std::endl;

        }

        config::printline(MvsTinfofs);
    }
    if(MvsTCase4==true)
    {
        // Case 4 - No dipolar, with thermal, with exchange
        config::printline(MvsTinfofs);
        MvsTinfofs << "Running MvsT case 4" << std::endl;
        fields::checkdipolar=false;
        fields::checkthermal=true;
        fields::checkexchange=true;
        tdp::exchstifffp=tdp::exchstiff0;
        MvsTopfs << std::setprecision(16);
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            spins::sx[i]=0.1;
            spins::sy[i]=0.0;
            spins::sz[i]=sqrt(1.0-spins::sx[i]*spins::sx[i]);
        }
        MvsTopfs << std::endl << std::endl;
        MvsTopfs << "#four" << std::endl;
        for(double T = MvsTminTemp ; T < MvsTmaxTemp ; T+=MvsTdTemp)
        {
            //set the uniform system temperature
            tdp::uniformtemp=T;
            util::RunningStat MC;
            MC.Clear();
            MvsTinfofs << std::endl;
            FIXOUT(MvsTinfofs,"Converging temperatue:" << T << std::endl << std::endl);
            double oldmean=0.0;
            #ifdef CUDA
            cullb::initGPU();
            #endif

            for(unsigned int t = 0 ; ;)
            {
                util::calcm();
                MC.Push(spins::modm);
                if(t>0 && t%1000==0)
                {
                    MvsTinfofs.width(15);MvsTinfofs << "| Mean = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Mean() << " | delta M = " << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << fabs(MC.Mean()-oldmean) << " [ " << convmean << " ] | Variance =" << std::showpos << std::fixed << std::setfill(' ') << std::setw(18) << MC.Variance() << " [ " << convvar << " ] " << std::endl;
                    if(((fabs(MC.Mean()-oldmean)) < convmean) && (MC.Variance()<convvar))
                    {
                        MvsTopfs << T << "\t" << spins::modm << std::endl;
                        break;
                    }
                    oldmean=MC.Mean();
                }
                LLB::integrate(t);
    //            std::cout << t << "\t" << spins::mx << "\t" << spins::my << "\t" << spins::mz << "\t" << spins::modm << std::endl;

            }
            MvsTinfofs << std::endl;

        }

        config::printline(MvsTinfofs);
    }

    MvsTopfs.close();
    if(MvsTopfs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close MvsT file");
    }

    if(boltz1d==true)
    {
        //=====================================TESTING THE BOLTZMANN DISTRIBUTION=======================================
        fields::checkthermal=true;
        fields::checkdipolar=false;
        fields::checkexchange=false;
        tdp::exchstifffp=tdp::exchstiffNULL;
        fields::anisfp=fields::anis0;
        tdp::uniformtemp=b1dtemp;
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            maths::sphereRandom(spins::sx[i],spins::sy[i],spins::sz[i]);
        }
        unsigned int bin[1001];
        for(unsigned int i = 0 ; i < 1001 ; i++){bin[i]=0;}
        unsigned long int counter=0;
        #ifdef CUDA
        cullb::initGPU();
        #endif

        for(unsigned int t = 0 ; t < b1dts ; )
        {
            LLB::integrate(t);


            if(t>10000)
            {
                for(unsigned int i = 0 ; i < geom::ss ; i++)
                {
                    double modm=sqrt(spins::sx[i]*spins::sx[i]+spins::sy[i]*spins::sy[i]+spins::sz[i]*spins::sz[i]);
                    bin[int(modm*1000+0.5)]++;
                    counter++;
                }
            }
        }
        bfofs << "#|m| \t p(|m|)" << std::endl;
        for(unsigned int i = 0 ; i < 1001; i++)
        {
            bfofs << double(i)/1000.0 << "\t" << bin[i]/(double(counter)) << std::endl;
        }
    }


}
