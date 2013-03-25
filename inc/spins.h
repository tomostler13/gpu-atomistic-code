// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 25 Mar 2013 17:25:14
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
    extern Array3D<double> Sznzp;
    extern Array3D<fftw_complex> Sqznzp;
    extern unsigned int nzpcplxdim;
    extern fftw_plan SzcfPF,SzcfPB;
    extern double normsize;
	extern unsigned int update;
    void initSpins(int argc,char *argv[]);
    void FFTForward();
    void eFFTForward();
    void calcRealSpaceCorrelationFunction(unsigned int);
    extern std::string rscfstr;
    extern std::ofstream rscfs;
}
#endif /*_SPINS_H_*/
