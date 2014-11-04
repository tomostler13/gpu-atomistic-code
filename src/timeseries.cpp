// File: timeseries.cpp
// Author: Tom Ostler
// Created: 03 Nov 2014
// Last-modified: 04 Nov 2014 17:13:39
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
    //This array stores the spin map in real space
    Array3D<fftw_complex> s3d;
    s3d.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
    s3d.IFill(0);
    fftw_plan ftspins;
    //plan an in-place transform
    fftw_set_timelimit(60);
    ftspins=fftw_plan_dft_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],s3d.ptr(),s3d.ptr(),FFTW_FORWARD,FFTW_PATIENT);
    SUCCESS(config::Info);
    num_samples=dsf::rts/(spins::update*dsf::dsfupdate);
    FIXOUT(config::Info,"Number of time samples:" << num_samples << std::endl);
    Array2D<fftw_complex> Cqt;
    Cqt.resize(dsf::nk,num_samples);



    unsigned int sample_counter=0;
    for(unsigned int t = 0 ; t < dsf::ets ; t++)
    {
        if(t%spins::update==0)
        {
            util::calc_mag();
            util::output_mag(magout,t);
        }
        llg::integrate(t);
    }
    std::ofstream
    for(unsigned int k = 0 ; k < dsf::nk ; k++)
    {
        int kvec[3]={dsf::kpoints(k,0),dsf::kpoints(k,1),dsf::kpoints(k,2)};
        std::stringstream sstr;
        sstr << "kx" << kvec[0] << "ky" << kvec[1] << "kz" << kvec[2] << ".dat";
        std::string str=sstr.str();
        std::ofstream kvout;
        kvout.open(str.c_str());
        bool fopen=false;
        if(!kvout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not open file for outputting PSD for a k-vector. Redirecting to cout");
        }
        else
        {
            fopen=true;
        }
        //output the sampling information
        kvout << "#numsamples\t" << num_samples << std::endl;
        kvout << "#samplestime\t" << static_cast<double>(dsf::dsfupdate*spins::update)*llg::dt << std::endl;
        kvout << "#kvec\t" << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << std::endl;
        kvout << "#dim\t" << geom::dim[0] << "\t" << geom::dim[1] << "\t" << geom::dim[2] << std::endl;
        kvout << "#kpointdim\t" << geom::dim[0]*geom::Nk[0] << "\t" << geom::dim[1]*geom::Nk[1] << "\t" << geom::dim[2]*geom::Nk[2] << std::endl;
            int negq=-static_cast<int>(i)+num_samples;
            double freq=(static_cast<double>(i)*2.0*M_PI)/(static_cast<double>(dsf::dsfupdate*spins::update*num_samples)*llg::dt);
            //calculate the bits of the one-sided PSD
            const double Hq = Cqt(k,i)[0]*Cqt(k,i)[0] + Cqt(k,i)[1]*Cqt(k,i)[1];
            const double Hmq = Cqt(k,negq)[0]*Cqt(k,negq)[0] + Cqt(k,negq)[1]*Cqt(k,negq)[1];
            if(fopen)
            {
                kvout << i << "\t" << freq << "\t" << Hq+Hmq << std::endl;
            }
            else//try to redirect cout
            {
                std::cout << i << "\t" << freq << "\t" << Hq+Hmq << std::endl;
            }
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

    for(unsigned int t = dsf::ets ; t < (dsf::rts+dsf::ets) ; t++)
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
            fftw_execute(ftspins);
            //loop over the k-points that we are interested in
            for(unsigned int i = 0 ; i < dsf::nk ; i++)
            {
                int kvec[3]={dsf::kpoints(i,0),dsf::kpoints(i,1),dsf::kpoints(i,2)};
                int relkvec[3]={0,0,0};
                for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                {
                    if(kvec[xyz]<0)
                    {
                        relkvec[xyz]=kvec[xyz]+geom::dim[xyz]*geom::Nk[xyz];
                    }
                }
                if(sample_counter<num_samples)
                {
                    Cqt(i,sample_counter)[0]=s3d(relkvec[0],relkvec[1],relkvec[2])[0];
                    Cqt(i,sample_counter)[1]=s3d(relkvec[0],relkvec[1],relkvec[2])[1];
                }
            }
            sample_counter++;

        }
        llg::integrate(t);
    }
}
