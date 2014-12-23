// File: timeseries.cpp
// Author: Tom Ostler
// Created: 03 Nov 2014
// Last-modified: 09 Dec 2014 20:36:36
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
void sim::timeseries(int argc,char *argv[])
{
    unsigned int num_samples=0;
    double et=0.0,rt=0.0;
    //ets = equilibration time
    //rts = run timesteps
    unsigned int ets=0,rts=0;
    if(sf::csf==false)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The flag sf:Calculate must be set for true so this simulation.");
    }
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Time series details***" << std::endl;
    //The k-vector information is read in the file sf.cpp and sf_glob.cpp
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

    libconfig::Setting &setting = config::cfg.lookup("timeseries");
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
    errstatus=setting.lookupValue("Temperature",llg::T);
    if(errstatus)
    {
        FIXOUT(config::Info,"Temperature:" << llg::T << " [K]" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read temperature for sf (sf:Temperature (double))");
    }
    bool oits=false;
    errstatus=setting.lookupValue("OutputIndividualTimeSeries",oits);
    if(errstatus)
    {
        FIXOUT(config::Info,"Outputting individual time series for each magnetic species?:" << config::isTF(oits) << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read whether you want to output the individual time series for each magnetic species (timseries:OutputIndividualTimeSeries (bool))");
    }

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


    //calculate the total number of samples
    num_samples=rts/(spins::update*sf::sfupdate);
    FIXOUT(config::Info,"Number of time samples:" << num_samples << std::endl);
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
    kvinfo << "#numsamples\t" << num_samples << std::endl;
    kvinfo << "#samplestime\t" << static_cast<double>(sf::sfupdate*spins::update)*llg::dt << std::endl;
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
            ikvinfo[i]<< "#numsamples\t" << num_samples << std::endl;
            ikvinfo[i]<< "#samplestime\t" << static_cast<double>(sf::sfupdate*spins::update)*llg::dt << std::endl;
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



    for(unsigned int t = 0 ; t < ets ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
        }
        llg::integrate(t);
    }

    unsigned int sample_counter=0;



    spins::mapout=true;
    for(unsigned int t = ets ; t < (rts+ets) ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(t);
        }
        if(t%(spins::update*sf::sfupdate)==0)
        {

                if(sample_counter<num_samples)
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

                    //output the info for each magnetic sublattice
                    is3d.IFill(0);
                    for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
                    {
                        //copy spin to 3D arrays and DO NOT apply unitary transforms
                        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            unsigned int xyz[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                            //get the magnetic species number
                            unsigned int ms=geom::lu(i,3);
                            if(ms==s)//then output store the info for that magetic species
                            {
                                is3d(s,xyz[0],xyz[1],xyz[2])[0]=spins::Sx[i];
                                is3d(s,xyz[0],xyz[1],xyz[2])[1]=spins::Sy[i];
                            }
                            else
                            {
                                is3d(s,xyz[0],xyz[1],xyz[2])[0]=0.0;
                                is3d(s,xyz[0],xyz[1],xyz[2])[1]=0.0;
                            }
                        }
                    }
                    //execute the plan many (one for each species)
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
                }
                sample_counter++;

        }
        llg::integrate(t);
    }
    if(sample_counter!=num_samples)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("The number of samples counted (not output) is not the same as num_samples that was calculated before hand.");
    }
    if(kvout.is_open())
    {
        kvout.close();
        if(kvout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close k-vector file.");
        }
    }
    if(kvinfo.is_open())
    {
        kvinfo.close();
        if(kvinfo.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close k-vector info file.");
        }
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
