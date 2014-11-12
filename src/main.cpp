// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 06 Nov 2014 11:54:05
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
#ifdef CUDA
#include "../inc/cuda.h"
#endif
int main(int argc,char *argv[])
{
    config::initConfig(argc,argv);
    //Initialize the lattice
    geom::initGeom(argc,argv);
    if(config::exchm==0)
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
    //Read in the exchange matrix
    exch::initExch(argc,argv);

    //Now we have all of the terms in our interaction matrix, fourier transform the result
    if(config::dipm==0 && config::exchm==0)
    {
        intmat::fftIntmat();
    }
    else if(config::dipm==0 && config::inc_dip==true)
    {
        intmat::fftDipIntmat();
    }
    //Initialise the field arrays
    fields::initFields(argc,argv);
    //Initialise the spin arrays
    spins::initSpins(argc,argv);
    sim::initSim(argc,argv);
    llg::initLLG(argc,argv);
    //Initialise the Dynamic structure factor calculation
    sf::initSF(argc,argv);

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
    else if(sim::sim_type=="quick")
    {

        llg::T=1e-27;
        int counter=0;
        time_t now = time(0);
        char *dtime=ctime(&now);
        std::cout << "#Start time:\t" << dtime << std::endl;
        for(unsigned int t = 0 ; t < 10000 ; t++)
        {
            if(t%spins::update==0)
            {
                if(counter%10==0)
                {
                    counter=0;
                }
                counter++;
                util::calc_mag();
                std::cout << static_cast<double>(t)*llg::dt << "\t";
                for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
                {
                   std::cout << spins::mag(s,0) << "\t" << spins::mag(s,1) << "\t" << spins::mag(s,2) << "\t";
                }
                std::cout << std::endl;
            }

            llg::integrate(t);
        }
        time_t nowe=time(0);
        char *dtimee=ctime(&nowe);
        std::cout << "#End time:\t" << dtimee << std::endl;
    }
    return(0);
}
