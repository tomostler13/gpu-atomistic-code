// File: timeseries.cpp
// Author: Tom Ostler
// Created: 03 Nov 2014
// Last-modified: 04 Nov 2014 20:37:25
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
#include "../inc/dsf.h"
void sim::timeseries(int argc,char *argv[])
{
    unsigned int num_samples=0;
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Timeseriesdetails***" << std::endl;
    //at the moment all of the information for the timeseries (to calculate the dsf)
    //is pre read in the files dsf.cpp dsf_glob.cpp. Eventually I would like to
    //integrate the DSF so that it can be run as part of any simulation.

    //set the simulation temperature
    llg::T=dsf::T;
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



    //calculate the total number of samples
    num_samples=dsf::rts/(spins::update*dsf::dsfupdate);
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


    for(unsigned int t = 0 ; t < dsf::ets ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(magout,t);
        }
        llg::integrate(t);
    }

    unsigned int sample_counter=0;


    //output the sampling information
    kvinfo << "#datafile\t" << "kvec.dat" << std::endl;
    kvinfo << "#numsamples\t" << num_samples << std::endl;
    kvinfo << "#samplestime\t" << static_cast<double>(dsf::dsfupdate*spins::update)*llg::dt << std::endl;
    kvinfo << "#dim\t" << geom::dim[0] << "\t" << geom::dim[1] << "\t" << geom::dim[2] << std::endl;
    kvinfo << "#kpointdim\t" << geom::dim[0]*geom::Nk[0] << "\t" << geom::dim[1]*geom::Nk[1] << "\t" << geom::dim[2]*geom::Nk[2] << std::endl;
    kvinfo << "#NumberKPoints\t" << dsf::nk << std::endl;

    Array2D<unsigned int> kplu;
    kplu.resize(dsf::nk,3);
    for(unsigned int k = 0 ; k < dsf::nk ; k++)
    {
        int kvec[3]={dsf::kpoints(k,0),dsf::kpoints(k,1),dsf::kpoints(k,2)};
        for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
        {
            if(kvec[xyz]<0)
            {
                kplu(k,xyz)=kvec[xyz]+geom::dim[xyz]*geom::Nk[xyz];
            }
        }
        kvinfo << "#kvector " << k << " = " << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << std::endl;;
    }

    for(unsigned int t = dsf::ets ; t < (dsf::rts+dsf::ets) ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(magout,t);
        }
        if(t%(spins::update*dsf::dsfupdate)==0)
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
                        s3d(xyz[0],xyz[1],xyz[2])[0]=spins::Sx[i]*dsf::uo(ms,0);
                        s3d(xyz[0],xyz[1],xyz[2])[1]=spins::Sy[i]*dsf::uo(ms,1);
                    }
                    fftw_execute(ftspins);
                    //output the time for completeness
                    kvout << static_cast<double>(t-dsf::ets)*llg::dt << "\t";
                    //loop over the k-points that we are interested in
                    for(unsigned int i = 0 ; i < dsf::nk ; i++)
                    {
                        kvout << s3d(kplu(i,0),kplu(i,1),kplu(i,2))[0] << "\t" << s3d(kplu(i,0),kplu(i,1),kplu(i,2))[1] << "\t";
                    }
                    kvout << std::endl;
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

}
