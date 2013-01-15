// File: main.cpp
// Author:Tom Ostler
// Last-modified: 09 Jan 2013 11:20:28
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
//------header files for globals-------
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/mat.h"
#include "../inc/intmat.h"
#include "../inc/spins.h"
#include "../inc/fields.h"
#include "../inc/maths.h"
#include "../inc/util.h"
#include "../inc/neigh.h"
#include "../inc/LLBCPU.h"
#include "../inc/LLB.h"
#include "../inc/tdp.h"
#include "../inc/sim.h"
#include "../inc/mf.h"
#ifdef CUDA
#include "../inc/cuda.h"
void cullb::cuinit(int argc,char *argv[]);
#endif /*CUDA*/
//-------------------------------------
#define FIXOUT(a,b) a.width(75);a << std::left << b;
int main(int argc,char *argv[])
{

	//setup the reading of the configuration file and
	//print out some basic stuff (compiler,date,git sha etc)
    config::initConfig(argc,argv);
	//lookup to geometry (system dimensions etc)
    geom::initGeom(argc,argv);
	//initialize the material properties (Tc,Ms...)
	mat::initMat(argc,argv);
    //initialize the mean field arrays
    mf::initmf(argc,argv);
	//initialize fields
	fields::initFields(argc,argv);

    tdp::inittdp(argc,argv);


	assert(fields::fi);
	if(fields::dfc=="fft")
	{
        FIXOUT(config::Info,"Using FFT for dipolar field calculations:" << "TRUE" << std::endl);
		//setup arrays for interaction matrix
		intmat::initIntmat(argc,argv);
		//fill the interaction matrix for the demag field calculations
		intmat::fillIntmat();
        //initialise the transforms
        intmat::fftIntmat();
	}
    //initialise the neighour list
    neigh::initNeigh(argc,argv);

    //initialise the utilities stuff
    util::initUtil(argc,argv);
	//initialise the spin arrays
	spins::initSpins(argc,argv);

    sim::initSim(argc,argv);
    #ifdef CUDA
    cullb::cuinit(argc,argv);
    #endif
    llb::initLLBCPU(argc,argv);

    if(sim::sim_type=="test temp dep prop")
    {
        sim::testtempdepprop(argc,argv);
    }
    else if(sim::sim_type=="test suite")
    {
        sim::testsuite(argc,argv);
    }
    else if(sim::sim_type=="quick")
    {
        tdp::uniformtemp=1;
        tdp::Tcalc();
        #ifdef CUDA
        cullb::initGPU();
        #endif /*CUDA*/
        int counter=0;
        for(unsigned int t=0 ; t < 50000000 ; )
        {
//            if((t%(util::update))==0 || t==0)
//            {
                util::outputSpinsVTU(t);
//            }

            util::calcm();
            LLB::integrate(t);
            std::cout << t << "\t" << spins::mx << "\t" << spins::my << "\t" << spins::mz << "\t" << spins::modm << std::endl;
        }
        //output final spin state
        util::outputSpinsVTU(-1);
    }
    else if(sim::sim_type=="nothing")
    {
        //do nothing
        FIXOUT(config::Info,"Doing nothing:" << std::flush);
        config::Info << "that was easy" << std::endl;
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Phew, I thought I was going to have to do some work then but I don't recognise the simulation type. Please don't try again.");
    }

    return(0);
}
