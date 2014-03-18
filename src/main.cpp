// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 18 Mar 2014 14:24:55
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
#include "../inc/mat.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/exch.h"
#include "../inc/anis.h"
#include "../inc/llgCPU.h"
#include "../inc/util.h"
#include "../inc/sim.h"
#include <omp.h>
#include "../inc/llg.h"
#ifdef CUDA
#include "../inc/cuda.h"
#endif
int main(int argc,char *argv[])
{
	config::initConfig(argc,argv);
	//Initialize the lattice
	geom::initGeom(argc,argv);
	//Read the material properties
	mat::initMat(argc,argv);
	//initialize the interaction matrices
	intmat::initIntmat(argc,argv);

	//Read in the exchange matrix
	exch::initExch(argc,argv);
	//Read in the anisotropy tensor
	anis::initAnis(argc,argv);
	//add the dipolar fields
	intmat::fillIntmat();


	//Now we have all of the terms in our interaction matrix, fourier transform the result
	intmat::fftIntmat();
	//Initialise the field arrays
	fields::initFields(argc,argv);
	//Initialise the spin arrays
	spins::initSpins(argc,argv);
	sim::initSim(argc,argv);
	llg::initLLG(argc,argv);
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
	else if(sim::sim_type=="quick")
	{

		llg::T=1.0e-27;
        int counter=0;
		for(unsigned int t = 0 ; t < 50000 ; t++)
		{
			llg::integrate(t);
			if(t%spins::update==0)
			{
                if(counter%10==0)
                {
                    util::outputSpinsVTU(t);
                    counter=0;
                }
                counter++;
				const double mx = util::reduceCPU(spins::Sx,geom::nspins);
				const double my = util::reduceCPU(spins::Sy,geom::nspins);
				const double mz = util::reduceCPU(spins::Sz,geom::nspins);
				std::cout << double(t)*llg::dt << "\t" << mx/double(geom::nspins) << "\t" << my/double(geom::nspins) << "\t" << mz/double(geom::nspins) << std::endl;
			}
		}
        util::outputSpinsVTU(-1);
	}
    return(0);
}
