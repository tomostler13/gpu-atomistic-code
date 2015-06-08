// File: laser_heating.cpp
// Author: Tom Ostler
// Created: 24 Nov 2014
// Last-modified: 08 Jun 2015 12:04:53
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
void CalcTe(double*,double*,double,unsigned int,double,double,double&,double&,double,double,double,double,double,double);
void sim::laser_heating(int argc,char *argv[])
{
    bool opsf=false;//OutPut Structure Factor info
    unsigned int ets = 0 , rts = 0 , num_pulses=0;
    double et = 0.0 , rt = 0.0 , gamma_e = 0.0 , Cl = 0.0 , G_ep = 0.0 , Pump_Time = 0.0 , ct = 0.0 , pumpfluence = 0 , initT = 0.0;
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
    if(setting.lookupValue("OuputIndividualTimeSeries",oits))
    {
        FIXOUT(config::Info,"Outputting individual ssf's for each species?:" << config::isTF(oits) << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read setting laser_heating.OutputIndividualTimeSeries (bool)");
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
    errstatus = setting.lookupValue("PumpFluence",pumpfluence);
    if(errstatus)
    {
        FIXOUT(config::Info,"Laser power:" << pumpfluence << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the pump fluence (laserheating:PumpFluence (double))");
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

    FIXOUT(config::Info,"Planning FFT for spin map:" << std::flush);
    //This array stores the spin map in real space
    Array3D<fftw_complex> s3d;
    //this is the same as s3d except that it is for calculating the structure factor for each magnetic species
    Array4D<fftw_complex> is3d;
    s3d.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
    s3d.IFill(0);
    fftw_plan ftspins;
    //plan an in-place transform
    fftw_set_timelimit(60);
    ftspins=fftw_plan_dft_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],s3d.ptr(),s3d.ptr(),FFTW_FORWARD,FFTW_PATIENT);
    SUCCESS(config::Info);
    fftw_plan iftspins;
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
    //we are declaring these on the heap because with clang 6.0 (at least)
    //the declaration on the stack is not allowed. I didn't manage to work
    //out why this is the case.
    std::ofstream *ikvinfo,*ikvout;

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
            util::output_mag(t);
            ttmout << static_cast<double>(t)*llg::dt << "\t" << initT << "\t" << initT << std::endl;
        }
        llg::integrate(t);
    }


    // Check the flag to output the spin map
    spins::mapout=true;
    double Te=initT,Tl=initT;
    for(unsigned int t = ets ; t < ets+rts ; t++)
    {
        unsigned int nt=t-ets;
        CalcTe(pulse_scale,pulse_delays,static_cast<double>(nt)*llg::dt,num_pulses,Pump_Time,pumpfluence,Te,Tl,G_ep,Cl,gamma_e,llg::dt,initT,ct);
        llg::T=Te;
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
            ttmout << static_cast<double>(t)*llg::dt << "\t" << Te << "\t" << Tl << std::endl;
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
        if(oits)
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
                        if(ms==s)
                        {
                            if(sf::qa[0]>1e-12)//quantization axis is x
                            {
                                is3d(s,xyz[0],xyz[1],xyz[2])[0]=spins::Sy[i]*sf::uo(ms,1);
                                is3d(s,xyz[0],xyz[1],xyz[2])[1]=spins::Sz[i]*sf::uo(ms,2);
                            }
                            else if(sf::qa[1]>1e-12)
                            {
                                is3d(s,xyz[0],xyz[1],xyz[2])[0]=spins::Sx[i]*sf::uo(ms,0);
                                is3d(s,xyz[0],xyz[1],xyz[2])[1]=spins::Sz[i]*sf::uo(ms,2);
                            }
                            else if(sf::qa[2]>1e-12)
                            {
                                is3d(s,xyz[0],xyz[1],xyz[2])[0]=spins::Sx[i]*sf::uo(ms,0);
                                is3d(s,xyz[0],xyz[1],xyz[2])[1]=spins::Sy[i]*sf::uo(ms,1);
                            }
                            else
                            {
                                is3d(s,xyz[0],xyz[1],xyz[2])[0]=0.0;
                                is3d(s,xyz[0],xyz[1],xyz[2])[1]=0.0;
                            }
                        }
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
