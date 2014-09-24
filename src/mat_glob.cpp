// File: mat_glob.cpp
// Author:Tom Ostler
// Created: 26 June 2014
// Last-modified: 23 Sep 2014 13:13:50
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
    Array<double> lambda;
    //Gyromagnetic ratio
    Array<double> gamma;
    //Bohr magneton
    double muB=9.27e-24;
    //Gyromagnetic ratio
    double gyro=1.76e11;
    //magnetic moment (muB)
    Array<double> mu;
    //thermal term prefactor
    Array<double> sigma;
}
