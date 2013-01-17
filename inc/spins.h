// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 17 Jan 2013 14:00:53
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../inc/arrays.h"
#ifndef _SPINS_H_
#define _SPINS_H_
namespace spins
{
    extern Array3D<fftw_complex> Skx,Sky,Skz;
    extern Array<double> Sx,Sy,Sz;
    void initSpins(int argc,char *argv[]);
    void FFTForward();
}
#endif /*_SPINS_H_*/
