// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 21 Jan 2013 16:01:30
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/intmat.h"
#include "../inc/mat.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/exch.h"
#include "../inc/anis.h"
#include "../inc/llg.h"
#include "../inc/util.h"
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
//    exch::initExch(argc,argv);
    //Read in the anisotropy tensor
    anis::initAnis(argc,argv);


//    intmat::fillIntmat();
    intmat::fftIntmat();

    //Initialise the field arrays
    fields::initFields(argc,argv);
    //Initialise the spin arrays
    spins::initSpins(argc,argv);
    llg::initLLG(argc,argv);
    //fields::bfdip();
    //fields::ftdip();
    llg::T=1.0e-27;
    for(unsigned int t = 0 ; t < 50000 ; t++)
    {
        llg::llgCPU(t);
        //const double mx = util::reduceCPU(sx,geom::nspins);
        const double mx = util::reduceCPU(spins::Sx,geom::nspins);
        const double my = util::reduceCPU(spins::Sy,geom::nspins);
        const double mz = util::reduceCPU(spins::Sz,geom::nspins);
        std::cout << mx/double(geom::nspins) << "\t" << my/double(geom::nspins) << "\t" << mz/double(geom::nspins) << std::endl;
    }
    return(0);
}
