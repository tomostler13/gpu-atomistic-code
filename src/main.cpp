// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 15 Jan 2013 19:55:10
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
int main(int argc,char *argv[])
{
    config::initConfig(argc,argv);
    //Initialize the lattice
    geom::initGeom(argc,argv);
    return(0);
}
