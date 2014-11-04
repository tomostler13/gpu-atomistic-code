// File: dsf.cpp
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 04 Nov 2014 12:49:17
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/dsf.h"
#include "../inc/unitcell.h"
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

//The aim of this part of the code is to take the user input
//and calculate the dynammic structure factor based on what
//the user selects. The actual calculations should not be
//specific to any one calculations and therefore have their
//own control features in the config file dsf:
namespace dsf
{
    void initDSF(int argc,char *argv[])
    {
        readDSFParam(argc,argv);
    }
}
