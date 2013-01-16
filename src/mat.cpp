// File: mat.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 16 Jan 2013 17:36:49
#include "../inc/mat.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/arrays.h"
#include <libconfig.h++>
#include <cassert>
#include <sstream>
#include <string>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace mat
{
    //bath coupling (no units)
    double lambda=0.0;
    //Gyromagnetic ratio
    double gamma=0.0;
    //Bohr magneton
    double muB=9.27e-24;
    //magnetic moment (muB)
    double mu=2.0;
    void initMat(int argc,char *argv[])
    {
        lambda=0.0;
        gamma=0.0;
        muB=9.27e-24;
        mu=2.0;
    }
}
