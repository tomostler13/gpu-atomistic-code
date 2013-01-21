// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 21 Jan 2013 15:16:27
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
    extern Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
    void initSpins(int argc,char *argv[]);
    void FFTForward();
    void eFFTForward();
}
#endif /*_SPINS_H_*/
