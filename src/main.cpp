// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 13 Sep 2016 11:35:04
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include <fftw3.h>
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/intmat.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/exch.h"
#include "../inc/anis.h"
#include "../inc/llgCPU.h"
#include "../inc/util.h"
#include "../inc/sim.h"
#include "../inc/llg.h"
#include "../inc/sf.h"
#include "../inc/rscf.h"
#include "../inc/defines.h"
#include "../inc/arrays.h"
#ifdef CUDA
#include "../inc/cuda.h"
#endif
int main(int argc,char *argv[])
{
    config::initConfig(argc,argv);
    //Initialize the lattice
    geom::initGeom(argc,argv);
    if(config::exchm==0)//then we are using the FFT method for the exchange calculation
    {
        //initialize the interaction matrices
        intmat::initIntmat(argc,argv);

        if(config::inc_dip==true && config::dipm==0)
        {
            //add the dipolar fields
            intmat::fillIntmat();
        }
    }
    else if(config::inc_dip==true && config::dipm==0)
    {
        intmat::initDipIntmat(argc,argv);
        // in this case the interaction matrix does
        // not need to be species dependent
        intmat::fillDipIntmat();
    }
    else if(config::exchm>98)
    {
        //If we are using the hybrid method we will always
        //use CSR for the CPU and DIA for the GPU
        intmat::initIntmat(argc,argv);
        if(config::inc_dip==true && config::dipm==0)
        {
            intmat::initDipIntmat(argc,argv);
            intmat::fillDipIntmat();
        }

    }
    //Read in the exchange matrix
    exch::initExch(argc,argv);
    if(exch::eaem)
    {
        exit(EXIT_SUCCESS);
    }
    //Now we have all of the terms in our interaction matrix, fourier transform the result
    if(config::dipm==0 && config::exchm==0)
    {
        intmat::fftIntmat();
    }
    else if(config::dipm==0 && config::inc_dip==true)
    {
        intmat::fftDipIntmat();
    }
    else if(config::exchm>98)
    {
        intmat::fftIntmat();
    }
    //Initialise the field arrays
    fields::initFields(argc,argv);
    //Initialise the spin arrays
    spins::initSpins(argc,argv);
    sim::initSim(argc,argv);
    llg::initLLG(argc,argv);
    //Initialise the Dynamic structure factor calculation
    sf::initSF(argc,argv);
    //initialise the real space correlation function calculations
    rscf::initRSCF(argc,argv);
    //unsigned int temp=0;
    //rscf::calcRSCF(temp);
#ifdef CUDA
    cullg::cuinit(argc,argv);
#else
    llgCPU::initLLG(argc,argv);
#endif

    if(sim::sim_type=="MvT")
    {
        sim::MvT(argc,argv);
    }
    else if(sim::sim_type=="suscep")
    {
        //        sim::suscep(argc,argv);
    }
    else if(sim::sim_type=="timeseries")
    {
        sim::timeseries(argc,argv);
    }
    else if(sim::sim_type=="laserheating")
    {
        sim::laser_heating(argc,argv);
    }
    else if(sim::sim_type=="fix_field")
    {
        sim::ramp_field(argc,argv);
    }
    else if(sim::sim_type=="thermal_hyst")
    {
        sim::thermal_hyst(argc,argv);
    }
    else if(sim::sim_type=="quick")
    {

        llg::T=0.1;
        int counter=0;
        time_t now = time(0);
        char *dtime=ctime(&now);
        std::cout << "#Start time:\t" << dtime << std::endl;
        for(unsigned int t = 0 ; t < 100000 ; t++)
        {
            if(t%spins::update==0)
            {
                if(counter%10==0)
                {
                    util::outputSpinsVTU(t);
                    counter=0;
                }
                counter++;
                util::calc_mag();
                util::output_mag(t);
                std::cout << static_cast<double>(t)*llg::dt << "\t";
                for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
                {
                   std::cout << spins::mag(s,0) << "\t" << spins::mag(s,1) << "\t" << spins::mag(s,2) << "\t";
                }
                std::cout << std::endl;
            }

            llg::integrate(t);
        }
//        util::outputSpinsVTU(-1);
        time_t nowe=time(0);
        char *dtimee=ctime(&nowe);
        unsigned int temp=-1;
        util::outputSpinsVTU(temp);
        std::cout << "Outputting fourier transform of spin map..." << std::flush;
//        Array5D<fftw_complex> lSk;
//        Array5D<double> lSr;
//        lSr.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
//        lSk.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
/*        int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
        int *inembed=n;
        int *onembed=n;
        int istride=1;
        int ostride=1;
        int odist=geom::zps;
        int idist=geom::zps;
        FIXOUT(config::Log,"Parameters entering into FFTW plan of spins (forward transform)" << std::endl);
        FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
        FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
        FIXOUT(config::Log,"How many (FFT's) = " << geom::ucm.GetNMS()*3 << std::endl);
        FIXOUT(config::Log,"Pointer of real space spins (Sr):" << lSr.ptr() << std::endl);
        FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
        FIXOUT(config::Log,"istride = " << istride << std::endl);
        FIXOUT(config::Log,"idist = " << idist << std::endl);
        FIXOUT(config::Log,"Pointer of reciprocal space spins (Sk):" << lSk.ptr() << std::endl);
        FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
        FIXOUT(config::Log,"ostride = " << ostride << std::endl);
        FIXOUT(config::Log,"odist = " << odist << std::endl);
        FIXOUT(config::Log,"Direction (sign) " << FFTW_FORWARD << std::endl);
        FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
        fftw_plan SP = fftw_plan_many_dft(3,n,geom::ucm.GetNMS()*3,lSr.ptr(),inembed,istride,idist,lSk.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);*/
        Array3D<fftw_complex> Skx,Sky,Skz;
        Array3D<double> Srx,Sry,Srz;
        Skx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
        Srx.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
        Sky.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
        Sry.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
        Skz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
        Srz.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
        fftw_plan SxP = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Srx.ptr(),Skx.ptr(),FFTW_ESTIMATE);
        fftw_plan SyP = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sry.ptr(),Sky.ptr(),FFTW_ESTIMATE);
        fftw_plan SzP = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Srz.ptr(),Skz.ptr(),FFTW_ESTIMATE);
        Srx.IFill(0);
        Sry.IFill(0);
        Srz.IFill(0);
        Skx.IFill(0);
        Sky.IFill(0);
        Skz.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            int ic[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Srx(ic[0],ic[1],ic[2])=spins::Sx[i];
            Sry(ic[0],ic[1],ic[2])=spins::Sy[i];
            Srz(ic[0],ic[1],ic[2])=spins::Sz[i];
        }
        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
        std::ofstream opkx("kx.dat");
        std::ofstream opky("ky.dat");
        std::ofstream opkz("kz.dat");
        for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                {
                    opkx << i << "\t" << j << "\t" << k << "\t" << Skx(i,j,k)[0] << "\t" << Skx(i,j,k)[1] << std::endl;
                    opky << i << "\t" << j << "\t" << k << "\t" << Sky(i,j,k)[0] << "\t" << Sky(i,j,k)[1] << std::endl;
                    opkz << i << "\t" << j << "\t" << k << "\t" << Skz(i,j,k)[0] << "\t" << Skz(i,j,k)[1] << std::endl;
                }
            }
        }
        opkx.close();
        opky.close();
        opkz.close();
        Srx.IFill(0);
        Sry.IFill(0);
        Srz.IFill(0);
        Skx.IFill(0);
        Sky.IFill(0);
        Skz.IFill(0);
        for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                if(geom::lu(i,4)==s1)
                {
                    int ic[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                    Srx(ic[0],ic[1],ic[2])=spins::Sx[i];
                    Sry(ic[0],ic[1],ic[2])=spins::Sy[i];
                    Srz(ic[0],ic[1],ic[2])=spins::Sz[i];
                }
            }
            fftw_execute(SxP);
            fftw_execute(SyP);
            fftw_execute(SzP);
            std::stringstream sstrx,sstry,sstrz;
            sstrx << "kx_" << s1 << ".dat";
            sstry << "ky_" << s1 << ".dat";
            sstrz << "kz_" << s1 << ".dat";
            std::string strx=sstrx.str(),stry=sstry.str(),strz=sstrz.str();
            opkx.open(strx.c_str());
            opky.open(stry.c_str());
            opkz.open(strz.c_str());
            for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                    {
                        opkx << i << "\t" << j << "\t" << k << "\t" << Skx(i,j,k)[0] << "\t" << Skx(i,j,k)[1] << std::endl;
                        opky << i << "\t" << j << "\t" << k << "\t" << Sky(i,j,k)[0] << "\t" << Sky(i,j,k)[1] << std::endl;
                        opkz << i << "\t" << j << "\t" << k << "\t" << Skz(i,j,k)[0] << "\t" << Skz(i,j,k)[1] << std::endl;
                    }
                }
            }
        }

        std::cout << "#End time:\t" << dtimee << std::endl;
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        std::cout << sim::sim_type << std::endl;
        error::errMessage("Simulation not recognised. Options are: \n-Mvt\n-suscep\n-timeseries\n-laserheating\n-fix_field\n-quick");
    }
    return(0);
}
