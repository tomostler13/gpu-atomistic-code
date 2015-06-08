// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 06 Jun 2015 17:54:57
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
    extern Array4D<fftw_complex> dipSk,dipSr,hSr,hSk;
    extern Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
    extern Array2D<double> mag;
    extern unsigned int update,mag_calc_method;
    extern bool output_mag,mapout;
    void initSpins(int argc,char *argv[]);
    void FFTForward();
    void hFFTForward();
    void eFFTForward();
    void heFFTForward();
    void dipFFTForward();
    void dipeFFTForward();
    void setSpinsConfig();
    void setSpinsRandom();
    void setSpinsChequer();
}
#endif /*_SPINS_H_*/
