// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 17 Jan 2013 14:51:07
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
int main(int argc,char *argv[])
{
    config::initConfig(argc,argv);
    //Initialize the lattice
    geom::initGeom(argc,argv);
    //Read the material properties
    mat::initMat(argc,argv);
    //initialize the interaction matrices
    intmat::initIntmat(argc,argv);
    //Initialise the field arrays
    fields::initFields(argc,argv);
    //Initialise the spin arrays
    spins::initSpins(argc,argv);
    //fields::bfdip();
    fields::ftdip();
    return(0);
}
