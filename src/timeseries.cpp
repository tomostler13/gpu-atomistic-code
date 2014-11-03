// File: timeseries.cpp
// Author: Tom Ostler
// Created: 03 Nov 2014
// Last-modified: 03 Nov 2014 14:35:39
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
    llg::T=dsf::T;
    std::ofstream magout("mag.dat");
    if(!magout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file mag.dat");
    }
    FIXOUT(config::Info,"Planning FFT for spin map:" << std::flush);
    Array3D<fftw_complex> s3d;
    s3d.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
    s3d.IFill(0);
    fftw_plan ftspins;
    //plan an in-place transform
    fftw_set_timelimit(60);
    ftspins=fftw_plan_dft_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],s3d.ptr(),s3d.ptr(),FFTW_FORWARD,FFTW_PATIENT);
    SUCCESS(config::Info);
    num_samples=dsf::rts/(spins::update*dsf::dsfupdate);
    FIXOUT(config::Info,"Number of samples:" << num_samples << std::endl);
    Array<fftw_complex> Cqt;
    Cqt.resize(num_samples);
    unsigned int sample_counter=0;
    for(unsigned int t = 0 ; t < dsf::ets ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(magout,t);
        }
        if(t%(spins::update*dsf::dsfupdate)==0)
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
            //execute the fftw
            fftw_execute(ftspins);
            //find the complex conjugate of the k-vectors we care about
            double cc[2]={s3d(0,0,6)[0],-s3d(0,0,6)[1]};
            //multiply the correlation function by it's cc and store
            double z_time_cz[2]={s3d(0,0,6)[0]*cc[0]-s3d(0,0,6)[1]*cc[1],s3d(0,0,6)[0]*cc[1]+s3d(0,0,6)[1]*cc[0]};

            Cqt(sample_counter)[0]=z_time_cz[0];
            Cqt(sample_counter)[1]=z_time_cz[1];

        }
        llg::integrate(t);
    }
    //At the end we want to perform the FFT in time
    fftw_plan fttime;
    fttime=fftw_plan_dft_1d(num_samples,Cqt.ptr(),Cqt.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(fttime);
    unsigned int halfsamp=num_samples/2;
    for(unsigned int i = 0 ; i < halfsamp ; i++)
    {
        int negfreq=-static_cast<int>(i)+num_samples;
        double freq=(static_cast<double>(i)*2.0*M_PI)/(static_cast<double>(dsf::dsfupdate*spins::update*num_samples)*llg::dt);
        std::cout << i << "\t" << freq << "\t" << Cqt(i)[0]*Cqt(i)[0]+Cqt(i)[1]*Cqt(i)[1] << "\t" << Cqt(negfreq)[0]*Cqt(negfreq)[0]+Cqt(negfreq)[1]*Cqt(negfreq)[1] << std::endl;
    }
}
