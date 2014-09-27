// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 26 Sep 2014 16:16:51
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
    extern Array5D<fftw_complex> Sk;
    extern Array5D<fftw_complex> Sr;
    extern Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
    extern Array2D<double> mag;
    extern unsigned int update,mag_calc_method;
    void initSpins(int argc,char *argv[]);
    void FFTForward();
    void eFFTForward();
}
#endif /*_SPINS_H_*/
