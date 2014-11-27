// File: laser_heating.cpp
// Author: Tom Ostler
// Created: 24 Nov 2014
// Last-modified: 27 Nov 2014 11:49:27
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
//function prototype
double pump_power(double);
void CalcTe(double*,double*,double,unsigned int,double,double,double&,double&,double,double,double,double,double,double);
void sim::laser_heating(int argc,char *argv[])
{
    bool opsf=false;//OutPut Structure Factor info
    unsigned int ets = 0 , rts = 0 , num_pulses=0;
    double et = 0.0 , rt = 0.0 , gamma_e = 0.0 , Cl = 0.0 , G_ep = 0.0 , Pump_Time = 0.0 , ct = 0.0 , T_el_max = 0 , initT = 0.0;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Laser Heating Simulation details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("laserheating");
    bool errstatus=setting.lookupValue("EquilibrationTime",et);
    if(errstatus)
    {
        FIXOUT(config::Info,"Equilibration time:" << et << " [s]" << std::endl);
        ets=static_cast<unsigned int>((et/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    }
    errstatus=setting.lookupValue("RunTime",rt);
    if(errstatus)
    {
        FIXOUT(config::Info,"Run time:" << rt << " [s]" << std::endl);
        rts=static_cast<unsigned int>((rt/llg::dt)+0.5);
        FIXOUT(config::Info,"Number of timesteps for run:" << rts << std::endl);
    }
    errstatus=setting.lookupValue("InitialTemperature",initT);
    if(errstatus)
    {
        FIXOUT(config::Info,"Initial Temperature:" << initT << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read temperature (laserheating:InitialTemperature (double))");
    }
    errstatus=setting.lookupValue("Cl",Cl);
    if(errstatus)
    {
        FIXOUT(config::Info,"Lattice specific heat:" << Cl << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read lattice specific heat (laserheating:Cl (double))");
    }
    errstatus = setting.lookupValue("PumpTime",Pump_Time);
    if(errstatus)
    {
        FIXOUT(config::Info,"Laser pump time (width):" << Pump_Time << " [s]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read laser pump time (laserheating:PumpTime (double)).");
    }
    errstatus = setting.lookupValue("G_ep",G_ep);
    if(errstatus)
    {
        FIXOUT(config::Info,"Electron-phonon coupling constant:" << G_ep << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the electron-phonon coupling constant (laserheating:G_ep (double)).");
    }

    errstatus = setting.lookupValue("Gamma_e",gamma_e);
    if(errstatus)
    {
        FIXOUT(config::Info,"gamma_e (Ce(Te)=gamma_e*Te):" << gamma_e << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the gamma_e (Ce(T)=gamma_e*Te) (laserheating:Gamma_e (double))");
    }
    errstatus = setting.lookupValue("CoolingTime",ct);
    if(errstatus)
    {
        FIXOUT(config::Info,"Cooling time:" << ct << " [s]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the cooling time (laserheating:CoolingTime (double))");
    }
    errstatus = setting.lookupValue("OutputStructureFactor",opsf);
    if(errstatus)
    {
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
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the state of the outputting of the structure factor (laserheating:OutputStructureFactor (bool))");
    }
    errstatus = setting.lookupValue("MaxElectronTemperature",T_el_max);
    if(errstatus)
    {
        FIXOUT(config::Info,"Input max electron temperature:" << T_el_max << std::endl);
        //convert to power
        T_el_max = pump_power(T_el_max);
        FIXOUT(config::Info,"Laser power:" << T_el_max << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the max electron temperature (laserheating:MaxElectronTemperature (double))");
    }
    errstatus = setting.lookupValue("NumPulses",num_pulses);
    if(errstatus)
    {
        FIXOUT(config::Info,"Number of pulses:" << num_pulses << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the number of pulses (laserheating:NumPulses (unsigned int))");
    }
    llg::T=initT;
    double pulse_delays[num_pulses];
    FIXOUT(config::Info,"Pulse delays:"<< "[ " << std::flush);
    for(unsigned int i = 0 ; i < num_pulses ; i++)
    {
        pulse_delays[i]=setting["PulseDelays"][i];
        if(i<(num_pulses-1))
        {
            config::Info << pulse_delays[i] << " , ";
        }
    }
    config::Info << pulse_delays[num_pulses-1] << " ]" << std::endl;
    double pulse_scale[num_pulses];
    FIXOUT(config::Info,"Pulse scaling factors:"<< "[ " << std::flush);
    for(unsigned int i = 0 ; i < num_pulses ; i++)
    {
        pulse_scale[i]=setting["PulseScale"][i];
        if(i<(num_pulses-1))
        {
            config::Info << pulse_scale[i] << " , ";
        }
    }
    config::Info << pulse_scale[num_pulses-1] << " ]" << std::endl;

    //open the magnetization file and do some error handling
    std::ofstream magout("mag.dat");
    if(!magout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file mag.dat");
    }
    FIXOUT(config::Info,"Planning FFT for spin map:" << std::flush);
    //This array stores the spin map in real space
    Array3D<fftw_complex> s3d;
    s3d.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
    s3d.IFill(0);
    fftw_plan ftspins;
    //plan an in-place transform
    fftw_set_timelimit(60);
    ftspins=fftw_plan_dft_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],s3d.ptr(),s3d.ptr(),FFTW_FORWARD,FFTW_PATIENT);
    SUCCESS(config::Info);
    std::ofstream kvout("kvec.dat"),kvinfo("kvecinfo.dat");
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
    Array2D<unsigned int> kplu;
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

    std::ofstream ttmout("ttm.dat");
    if(!ttmout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open ttm file (ttm.dat)");
    }

    for(unsigned int t = 0 ; t < ets ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(magout,t);
            ttmout << static_cast<double>(t)*llg::dt << "\t" << initT << "\t" << initT << std::endl;
        }
        llg::integrate(t);
    }

    double Te=initT,Tl=initT;
    for(unsigned int t = ets ; t < ets+rts ; t++)
    {
        unsigned int nt=t-ets;
        CalcTe(pulse_scale,pulse_delays,static_cast<double>(t)*llg::dt,num_pulses,Pump_Time,T_el_max,Te,Tl,G_ep,Cl,gamma_e,llg::dt,initT,ct);
        llg::T=Te;
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(magout,t);
            ttmout << static_cast<double>(t)*llg::dt << "\t" << initT << "\t" << initT << std::endl;
        }
        if(t%(spins::update*sf::sfupdate)==0)
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
        llg::integrate(t);
    }
    magout.close();
    if(magout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close the magnetization file (mag.dat)");
    }
    ttmout.close();
    if(ttmout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close the ttm file (ttm.dat)");
    }
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
}

double pump_power(double x)//we do not pass by reference here deliberately because T_el_max set
{
    return((-5.9408+0.00722*x+4.11859e-5*(x*x)-3.66025e-10*x*x*x)*1.152E20);
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
