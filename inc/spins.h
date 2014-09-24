// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 23 Jan 2013 10:44:46
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
    extern Array3D<double> Srx,Sry,Srz;
    extern Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
	extern unsigned int update;
    void initSpins(int argc,char *argv[]);
    void FFTForward();
    void eFFTForward();
}
#endif /*_SPINS_H_*/
