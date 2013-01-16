// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 16 Jan 2013 17:48:55
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
int main(int argc,char *argv[])
{
    config::initConfig(argc,argv);
    //Initialize the lattice
    geom::initGeom(argc,argv);
    //Read the material properties
    mat::initMat(argc,argv);
    //initialize the interaction matrices
    intmat::initIntmat(argc,argv);
    //Fill the interaction matrix with the demag term
    intmat::fillIntmat();
    return(0);
}
