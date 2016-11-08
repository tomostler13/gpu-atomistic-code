// File: laser_heating.cpp
// Author: Tom Ostler
// Created: 24 Nov 2014
// Last-modified: 18 Oct 2016 13:18:23
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>
#include "../inc/error.h"
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
#include "../inc/sf.h"
#include "../inc/rscf.h"
//function prototype
void CalcTe(double*,double*,double,unsigned int,double,double,double&,double&,double,double,double,double,double,double);
void sim::laser_heating()
{
    bool opsf=false;//OutPut Structure Factor info
    unsigned int ets = 0 , rts = 0 , num_pulses=0;
    double et = 0.0 , rt = 0.0 , gamma_e = 0.0 , Cl = 0.0 , G_ep = 0.0 , Pump_Time = 0.0 , ct = 0.0 , pumpfluence = 0 , initT = 0.0;
    unsigned int nfv=0;
    Array2D<double> field_val;
    Array<double> field_times;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Laser Heating Simulation details***" << std::endl;

    if(!config::cfg.exists("laserheating"))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting \"laserheating\" does not exist. Check your config file.");
    }
    libconfig::Setting &setting = config::cfg.lookup("laserheating");
    if(setting.lookupValue("EquilibrationTime",et))
    {
        FIXOUT(config::Info,"Equilibration time:" << et << " [s]" << std::endl);
        ets=static_cast<unsigned int>((et/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Equilibration time could not be read. Using default value of 50ps.");
        et=50e-12;
        FIXOUT(config::Info,"Equilibration time (default):" << et << " [s]" << std::endl);
        ets=static_cast<unsigned int>((et/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps (default):" << ets << std::endl);
    }

    if(setting.lookupValue("RunTime",rt))
    {
        FIXOUT(config::Info,"Run time:" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run:" << rts << std::endl);
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Run time could not be read. Using default value of 50ps.");
        rt=50e-12;
        FIXOUT(config::Info,"Run time (default):" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run (default):" << rts << std::endl);
    }
    if(setting.lookupValue("InitialTemperature",initT))
    {
        if(initT<1e-5)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Because the two temperature model has a 1/(Te*gamma_e) in it, if you set the value of the initial temperature to zero the temperature will blow up. Set the initial temperature (laserheating:InitialTemperature (double)) greater than 1e-5");
        }
        FIXOUT(config::Info,"Initial Temperature:" << initT << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read temperature (laserheating:InitialTemperature (double))");
    }
    bool oits=false;
    if(setting.lookupValue("OutputIndividualTimeSeries",oits))
    {
        FIXOUT(config::Info,"Outputting individual ssf's for each species?:" << config::isTF(oits) << std::endl);
        if(geom::Nkset==false && oits==true)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Can't presently calculate the ssf without specifying mesh.");
        }
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Outputting of invididual time series not set. Defaulting to false.");
        oits=false;
    }
    if(setting.lookupValue("Cl",Cl))
    {
        FIXOUT(config::Info,"Lattice specific heat:" << Cl << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read lattice specific heat (laserheating:Cl (double))");
    }
    if(setting.lookupValue("PumpTime",Pump_Time))
    {
        FIXOUT(config::Info,"Laser pump time (width):" << Pump_Time << " [s]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read laser pump time (laserheating:PumpTime (double)).");
    }
    if(setting.lookupValue("G_ep",G_ep))
    {
        FIXOUT(config::Info,"Electron-phonon coupling constant:" << G_ep << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the electron-phonon coupling constant (laserheating:G_ep (double)).");
    }

    if(setting.lookupValue("Gamma_e",gamma_e))
    {
        FIXOUT(config::Info,"gamma_e (Ce(Te)=gamma_e*Te):" << gamma_e << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the gamma_e (Ce(T)=gamma_e*Te) (laserheating:Gamma_e (double))");
    }
    if(setting.lookupValue("CoolingTime",ct))
    {
        FIXOUT(config::Info,"Cooling time:" << ct << " [s]" << std::endl);
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the cooling time (laserheating:CoolingTime (double)). Defaulting to a long time 3000e-9.");
        ct=3000e-9;
    }
    if(setting.lookupValue("OutputStructureFactor",opsf))
    {
        if(geom::Nkset==false && opsf==true)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You cannot output the structure factor without setting mesh points Nm");
        }
        FIXOUT(config::Info,"Output the information for the structure factor?" << config::isTF(opsf) << std::endl);
        if(opsf)
        {
            if(sf::csf==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("The flag laserheating:OutputStructureFactor was set to true and you have not set the structure factor and k-vectors up (sf:Calculate)");
            }
        }
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not read the state of the outputting of the structure factor (laserheating:OutputStructureFactor (bool)). Defaulting to false.");
    }

    if(setting.lookupValue("PumpFluence",pumpfluence))
    {
        FIXOUT(config::Info,"Laser power:" << pumpfluence << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the pump fluence (laserheating:PumpFluence (double))");
    }
    if(setting.lookupValue("NumPulses",num_pulses))
    {
        FIXOUT(config::Info,"Number of pulses:" << num_pulses << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the number of pulses (laserheating:NumPulses (unsigned int))");
    }
    llg::T=initT;
    double *pulse_delays=NULL;
    try
    {
        pulse_delays = new double[num_pulses];
    }
    catch(...)
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not malloc pulse_delays array");
    }
    FIXOUT(config::Info,"Pulse delays:"<< "[ " << std::flush);
    for(unsigned int i = 0 ; i < num_pulses ; i++)
    {
        try
        {
            pulse_delays[i]=setting["PulseDelays"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " must be set.";
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
        if(i<(num_pulses-1))
        {
            config::Info << pulse_delays[i] << " , ";
        }
    }
    config::Info << pulse_delays[num_pulses-1] << " ]" << std::endl;
    double *pulse_scale=new double [num_pulses];
    bool scale_pulses=false;
    if(setting.lookupValue("ScalePulses",scale_pulses))
    {
        FIXOUT(config::Info,"Scale pulses?:" << config::isTF(scale_pulses) << std::endl);
    }
    else
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not read if you wanted to scale the pulses, defaulting to false.");
    }
    if(scale_pulses)
    {
        FIXOUT(config::Info,"Pulse scaling factors:"<< "[ " << std::flush);
        for(unsigned int i = 0 ; i < num_pulses ; i++)
        {
            try
            {
                pulse_scale[i]=setting["PulseScale"][i];
            }
            catch(const libconfig::SettingNotFoundException &snf)
            {
                error::errPreamble(__FILE__,__LINE__);
                std::stringstream errsstr;
                errsstr << "Setting not found exception caught. Setting " << snf.getPath() << "\nYou have requested to scale the pulses but not provided a setting.";
                std::string errstr=errsstr.str();
                error::errMessage(errstr);
            }
            if(i<(num_pulses-1))
            {
                config::Info << pulse_scale[i] << " , ";
            }
        }
        config::Info << pulse_scale[num_pulses-1] << " ]" << std::endl;
    }
    else
    {
        for(unsigned int i = 0 ; i < num_pulses ; i++)
        {
            pulse_scale[i]=1.0;
        }
    }
    bool outVTU=false;
    double VTUstart=0.0;
    unsigned int VTUupdate=0;
    if(!setting.lookupValue("OutputSpinsVTU",outVTU))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not read whether you want to output the spin configurations to VTU files. laserheating:OutputSpinsVTU (bool), defaulting to false.");
        outVTU=false;
    }
    FIXOUT(config::Info,"Output spin config to VTU?:" << config::isTF(outVTU) << std::endl);
    if(outVTU)
    {
        //then read the parameters that it requires
        if(setting.lookupValue("VTUUpdate",VTUupdate))
        {
            FIXOUT(config::Info,"VTU Update:" << VTUupdate << " [demag]" << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the VTU update frequency. laserheating:VTUUpdate (int)");
        }
        if(setting.lookupValue("VTUOutputStart",VTUstart))
        {
            FIXOUT(config::Info,"Start of output of VTU files:" << VTUstart << " [s]" << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the start time of the VTU file output. laserheating:VTUOutputStart (double)");
        }
    }
    bool giveField=false;
    if(!setting.lookupValue("VaryField",giveField))
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not read if you want to vary the field (laser_heating.VaryField (bool), defaulting to false.");
        giveField=false;
    }
    FIXOUT(config::Info,"Vary field?:" << config::isTF(giveField) << std::endl);
    if(giveField)
    {
        if(setting.lookupValue("NumFieldValues",nfv))
        {
            FIXOUT(config::Info,"Number of field values:" << nfv << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read the number of field values. laserheating.NumFieldValues (int)");
        }
        field_val.resize(nfv,3);
        field_times.resize(nfv+1);
        field_val.IFill(0);field_times.IFill(0);
        for(unsigned int i = 0 ; i < nfv ; i++)
        {
            try
            {
                field_times[i+1]=setting["FieldTimes"][i];
            }
            catch(const libconfig::SettingNotFoundException &snf)
            {
                error::errPreamble(__FILE__,__LINE__);
                std::stringstream errsstr;
                errsstr << "Setting not found exception caught. Setting " << snf.getPath();
                std::string errstr=errsstr.str();
                error::errMessage(errstr);
            }
            std::stringstream sstr;
            sstr << "Field info -> value " << i+1 << " UP TO time " << field_times[i+1];
            std::string str=sstr.str();
            std::stringstream fsstr;
            fsstr << "FieldVal" << i+1;
            std::string fstr=fsstr.str();
            for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
            {

                try
                {
                    field_val(i,xyz)=setting[fstr.c_str()][xyz];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Setting " << snf.getPath();
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
            }
            FIXOUTVEC(config::Info,str,field_val(i,0),field_val(i,1),field_val(i,2));
        }
        config::Info << std::endl;
    }

    //declarations for the correlation function (may not be malloced)
    //This array stores the spin map in real space
    Array3D<fftw_complex> s3d;
    //this is the same as s3d except that it is for calculating the structure factor for each magnetic species
    Array4D<fftw_complex> is3d;
    fftw_plan ftspins,iftspins;
    std::ofstream kvout,kvinfo;
    Array2D<unsigned int> kplu;
    //we are declaring these on the heap because with clang 6.0 (at least)
    //the declaration on the stack is not allowed. I didn't manage to work
    //out why this is the case.
    std::ofstream *ikvinfo,*ikvout;

    if(sf::csf)
    {
        FIXOUT(config::Info,"Planning FFT for spin map:" << std::flush);
        s3d.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
        s3d.IFill(0);
        kvout.open("kvec.dat");
        if(!kvout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open kvec.dat");
        }
        kvinfo.open("kvecinfo.dat");
        if(!kvinfo.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open kvecinfo.dat");
        }
        //plan an in-place transform
        fftw_set_timelimit(300);
        ftspins=fftw_plan_dft_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],s3d.ptr(),s3d.ptr(),FFTW_FORWARD,FFTW_PATIENT);
        SUCCESS(config::Info);
        if(oits)
        {
            FIXOUT(config::Info,"Planning FFT for individual species:" << std::flush);
            is3d.resize(geom::ucm.GetNMS(),geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
            is3d.IFill(0);
            int rank=3;
            const int n[3] = {geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]};
            int howmany=geom::ucm.GetNMS();
            int idist=n[0]*n[1]*n[2],odist=n[0]*n[1]*n[2];
            int istride=1,ostride=1;
            const int *inembed=n,*onembed=n;
            iftspins=fftw_plan_many_dft(rank,n,howmany,is3d.ptr(),inembed,istride,idist,is3d.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_PATIENT);
            SUCCESS(config::Info);
        }
        if(!kvout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for outputting PSD for a k-vector.");
        }
        if(!kvinfo.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for outputting information on S(q,t).");
        }
        //output the sampling information
        kvinfo << "#datafile\t" << "kvec.dat" << std::endl;
        kvinfo << "#dim\t" << geom::dim[0] << "\t" << geom::dim[1] << "\t" << geom::dim[2] << std::endl;
        kvinfo << "#kpointdim\t" << geom::dim[0]*geom::Nk[0] << "\t" << geom::dim[1]*geom::Nk[1] << "\t" << geom::dim[2]*geom::Nk[2] << std::endl;
        kvinfo << "#NumberKPoints\t" << sf::nk << std::endl;

        kplu.resize(sf::nk,3);
        for(unsigned int k = 0 ; k < sf::nk ; k++)
        {
            int kvec[3]={static_cast<int>(sf::kpoints(k,0)),static_cast<int>(sf::kpoints(k,1)),static_cast<int>(sf::kpoints(k,2))};
            for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
            {
                if(kvec[xyz]<0)
                {
                    kplu(k,xyz)=kvec[xyz]+geom::dim[xyz]*geom::Nk[xyz];
                }
                else
                {
                    kplu(k,xyz)=kvec[xyz];
                }
            }
            kvinfo << "#kvector " << k << " = " << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << std::endl;;
            kvinfo << "#lookupkvector " << k << " = " << kplu(k,0) << "\t" << kplu(k,1) << "\t" << kplu(k,2) << std::endl;
        }
        //if we are outputting each (magnetic) sublattice time series then we need to open geom::ucm.GetNMS() files
        if(oits)
        {
            try
            {
                ikvinfo = new std::ofstream [geom::ucm.GetNMS()];
                ikvout = new std::ofstream [geom::ucm.GetNMS()];
            }
            catch(...)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not create the array std::ofstream instances on the heap (one for each species).");
            }

            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                std::stringstream sstr;
                sstr << "kvecinfo_spec_" << i << ".dat";
                std::string str=sstr.str();
                ikvinfo[i].open(str.c_str());
                if(!ikvinfo[i].is_open())
                {
                    sstr.str("");
                    sstr << "Could not open file for outputting info for k-vectors for magnetic species " << i;
                    str=sstr.str();
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage(str);
                }

                sstr.str("");
                sstr << "kvec_spec_" << i << ".dat";
                str=sstr.str();
                ikvout[i].open(str.c_str());
                if(!ikvout[i].is_open())
                {
                    sstr.str("");
                    sstr << "File kvec_spec_" << i << ".dat could not be opened.";
                    str=sstr.str();
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage(str);
                }
                //output the sampling information
                ikvinfo[i]<< "#datafile\t" << str << std::endl;
                ikvinfo[i]<< "#dim\t" << geom::dim[0] << "\t" << geom::dim[1] << "\t" << geom::dim[2] << std::endl;
                ikvinfo[i]<< "#kpointdim\t" << geom::dim[0]*geom::Nk[0] << "\t" << geom::dim[1]*geom::Nk[1] << "\t" << geom::dim[2]*geom::Nk[2] << std::endl;
                ikvinfo[i]<< "#NumberKPoints\t" << sf::nk << std::endl;
                for(unsigned int k = 0 ; k < sf::nk ; k++)
                {
                    int kvec[3]={static_cast<int>(sf::kpoints(k,0)),static_cast<int>(sf::kpoints(k,1)),static_cast<int>(sf::kpoints(k,2))};
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        if(kvec[xyz]<0)
                        {
                            kplu(k,xyz)=kvec[xyz]+geom::dim[xyz]*geom::Nk[xyz];
                        }
                        else
                        {
                            kplu(k,xyz)=kvec[xyz];
                        }
                    }
                    ikvinfo[i] << "#kvector " << k << " = " << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << std::endl;;
                    ikvinfo[i] << "#lookupkvector " << k << " = " << kplu(k,0) << "\t" << kplu(k,1) << "\t" << kplu(k,2) << std::endl;
                }
            }
        }
    }

    std::ofstream ttmout("ttm.dat");
    if(!ttmout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open ttm file (ttm.dat)");
    }
    else
    {
        ttmout << "#Time [s] - Te [K] - Tl [K] - Ts_1 - Ts_2.... - Ts_Nspec" << std::endl;
    }
    //counter for outputting of VTU files
    int VTUcount=0;
    //which field value are we on?
    int fv=0;
    //original field
    const double of[3]={llg::applied[0],llg::applied[1],llg::applied[2]};
    for(unsigned int t = 0 ; t < ets ; t++)
    {
        const double realtime=static_cast<double>(t)*llg::dt;
        if(giveField)
        {
            if(realtime >= field_times[fv] && realtime <= field_times[fv+1])
            {
                llg::applied[0]=of[0]+field_val(fv,0);
                llg::applied[1]=of[1]+field_val(fv,1);
                llg::applied[2]=of[2]+field_val(fv,2);
                fv++;
            }
        }
        else
        {
            llg::applied[0]=of[0];
            llg::applied[1]=of[1];
            llg::applied[2]=of[2];
        }
        //std::cout << realtime << "\t" << llg::applied[0] << "\t" << llg::applied[1] << "\t" << llg::applied[2] << std::endl;
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
            //This should really have an option to output less regularly
            rscf::calcRSCF(t);
        }
        if(t%spins::update==0)
        {
            if(outVTU)
            {
                if(VTUcount==static_cast<int>(VTUupdate))
                {
                    if(static_cast<double>(t)*llg::dt >= VTUstart)
                    {
                        util::outputSpinsVTU(t);
                    }
                    //restart the counter
                    VTUcount=0;
                }
                else
                {
                    VTUcount++;
                }
            }
        }
        llg::integrate(t);
        if(t%spins::update==0)
        {
            util::calc_Ts();
            ttmout << static_cast<double>(t)*llg::dt << "\t" << initT << "\t" << initT;
            for(unsigned int spec = 0 ; spec < geom::ucm.GetNMS() ; spec++)
            {
//                ttmout << "\t" << llg::cps(spec) << "\t" << llg::dps(spec) << "\t" << llg::Ts[spec];
                ttmout << "\t" << llg::Ts[spec];
            }
            ttmout << std::endl;
        }
    }


    // Check the flag to output the spin map
    spins::mapout=true;
    double Te=initT,Tl=initT;
    for(unsigned int t = ets ; t < ets+rts ; t++)
    {
        const double realtime=static_cast<double>(t)*llg::dt;
        if(giveField)
        {
            if(realtime >= field_times[fv] && realtime <= field_times[fv+1])
            {
                llg::applied[0]=of[0]+field_val(fv,0);
                llg::applied[1]=of[1]+field_val(fv,1);
                llg::applied[2]=of[2]+field_val(fv,2);
                fv++;
            }
        }
        else
        {
            llg::applied[0]=of[0];
            llg::applied[1]=of[1];
            llg::applied[2]=of[2];
        }
        //std::cout << realtime << "\t" << llg::applied[0] << "\t" << llg::applied[1] << "\t" << llg::applied[2] << std::endl;
        unsigned int nt=t-ets;
        CalcTe(pulse_scale,pulse_delays,static_cast<double>(nt)*llg::dt,num_pulses,Pump_Time,pumpfluence,Te,Tl,G_ep,Cl,gamma_e,llg::dt,initT,ct);
        llg::T=Te;
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
            rscf::calcRSCF(t);
        }
        if(t%spins::update==0)
        {
            if(outVTU)
                if(VTUcount==static_cast<int>(VTUupdate))
                {
                    if(static_cast<double>(t)*llg::dt >= VTUstart)
                    {
                        util::outputSpinsVTU(t);
                    }
                    //restart the counter
                    VTUcount=0;
                }
                else
                {
                    VTUcount++;
                }
        }
        if(opsf==true)
        {
            if(t%(spins::update*sf::sfupdate)==0 && opsf==true)
            {
                //zero the 3d spin array
                s3d.IFill(0);
                //copy spin to 3D arrays and apply unitary transforms
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    unsigned int xyz[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                    //get the magnetic species number
                    unsigned int ms=geom::lu(i,3);
                    s3d(xyz[0],xyz[1],xyz[2])[0]=spins::Sx[i]*sf::uo(ms,0);
                    s3d(xyz[0],xyz[1],xyz[2])[1]=spins::Sy[i]*sf::uo(ms,1);
                }
                fftw_execute(ftspins);
                //output the time for completeness
                kvout << static_cast<double>(t-ets)*llg::dt << "\t";
                //loop over the k-points that we are interested in
                for(unsigned int i = 0 ; i < sf::nk ; i++)
                {
                    kvout << s3d(kplu(i,0),kplu(i,1),kplu(i,2))[0] << "\t" << s3d(kplu(i,0),kplu(i,1),kplu(i,2))[1] << "\t";
                }
                kvout << std::endl;
            }

        }
        if(oits && opsf)
        {
            if(t%(spins::update*sf::sfupdate)==0)
            {
                is3d.IFill(0);
                for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
                {
                    for(unsigned int i = 0 ; i < geom::nspins; i++)
                    {
                        unsigned int xyz[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                        unsigned int ms=geom::lu(i,3);
                        is3d(s,xyz[0],xyz[1],xyz[2])[0]=spins::Sx[i]*sf::uo(ms,0);
                        is3d(s,xyz[0],xyz[1],xyz[2])[1]=spins::Sy[i]*sf::uo(ms,1);
                    }//spin loop
                }//species loop
                fftw_execute(iftspins);
                for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
                {
                    //output the time for completeness
                    ikvout[s] << static_cast<double>(t-ets)*llg::dt << "\t";
                    //loop over the k-points that we are interested in
                    for(unsigned int i = 0 ; i < sf::nk ; i++)
                    {
                        ikvout[s] << is3d(s,kplu(i,0),kplu(i,1),kplu(i,2))[0] << "\t" << is3d(s,kplu(i,0),kplu(i,1),kplu(i,2))[1] << "\t";
                    }
                    ikvout[s] << std::endl;
                }
            }//
        }//update time
        llg::integrate(t);
        if(t%spins::update==0)
        {
            util::calc_Ts();
            ttmout << static_cast<double>(t)*llg::dt << "\t" << Te << "\t" << Tl;
            for(unsigned int spec = 0 ; spec < geom::ucm.GetNMS() ; spec++)
            {
//                ttmout << "\t" << llg::cps(spec) << "\t" << llg::dps(spec) << "\t" << llg::Ts[spec];
                ttmout << "\t" << llg::Ts[spec];
            }
            ttmout << std::endl;
        }
    }
    ttmout.close();
    if(ttmout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close the ttm file (ttm.dat)");
    }
    if(sf::csf)
    {
        kvinfo.close();
        if(kvinfo.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close kvinfo.dat.");
        }
        kvout.close();
        if(kvout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close kvec.dat.");
        }
        if(oits)
        {
            for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
            {
                if(ikvout[s].is_open())
                {
                    ikvout[s].close();
                    if(ikvout[s].is_open())
                    {
                        std::stringstream sstr;
                        sstr << "Could not close the data file for species "  << s;
                        std::string str=sstr.str();
                        error::errPreamble(__FILE__,__LINE__);
                        error::errWarning(str);
                    }
                }
                if(ikvinfo[s].is_open())
                {
                    ikvinfo[s].close();
                    if(ikvinfo[s].is_open())
                    {
                        std::stringstream sstr;
                        sstr << "Could not close the info file for species "  << s;
                        std::string str=sstr.str();
                        error::errPreamble(__FILE__,__LINE__);
                        error::errWarning(str);
                    }
                }
            }
            try
            {
                delete [] ikvinfo;
                ikvinfo=NULL;
                delete [] ikvout;
                ikvout=NULL;
            }
            catch(...)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errWarning("Could not delete at least one of the species dependent info or data k-vec files.");
            }
        }
    }
    try
    {
        delete [] pulse_delays;
    }
    catch(...)
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not deallocate pulse_delays array.");
    }
    try
    {
        delete [] pulse_scale;
    }
    catch(...)
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not deallocate pulse_scale array.");
    }
}

void CalcTe(double *scale,double *delay,double time,unsigned int numpulses,double tau,double pumppower,double &Te,double &Tl,double G,double Cl,double gamma_e,double dt_real,double init_temp,double cool)
{
    double Ptot=0.0;
    //first calculate the power as a sum over the exponentials
    for(unsigned int i = 0 ; i < numpulses ; i++)
    {
        Ptot+=pumppower*scale[i]*exp(-((time-delay[i])/tau)*((time-delay[i])/tau));
    }

    Te = (-G*(Te-Tl)+Ptot)*dt_real/(gamma_e*Te) + Te;
    Tl = ( G*(Te-Tl)-Cl*(Tl-init_temp)/cool     )*dt_real/Cl + Tl;
}
