// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 08 Oct 2014 14:17:45
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
    extern Array4D<fftw_complex> dipSk,dipSr;
    extern Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
    extern Array2D<double> mag;
    extern unsigned int update,mag_calc_method;
    void initSpins(int argc,char *argv[]);
    void FFTForward();
    void eFFTForward();
    void dipFFTForward();
    void dipeFFTForward();
}
#endif /*_SPINS_H_*/
