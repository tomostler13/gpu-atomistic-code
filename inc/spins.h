// File: spins.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 19 Aug 2022 09:59:39
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
    extern Array3D<double> Srx,Sry,Srz;
    extern Array3D<fftw_complex> Skx,Sky,Skz;
    extern Array3D<fftw_complex> Ckxx,Ckxy,Ckxz,Ckyx,Ckyy,Ckyz,Ckzx,Ckzy,Ckzz,Cksds;
    extern Array3D<double> Crxx,Crxy,Crxz,Cryx,Cryy,Cryz,Crzx,Crzy,Crzz,Crsds;
    extern Array2D<double> mag;
    extern unsigned int update,mag_calc_method;
    extern bool output_mag,mapout;
    void initSpins();
    void FFTForward();
    void hFFTForward();
    void eFFTForward();
    void heFFTForward();
    void dipFFTForward();
    void dipeFFTForward();
    void setSpinsConfig();
    void setSpinsRandom();
    void setSpinsChequerX();
    void setSpinsChequerY();
    void setSpinsChequerZ();
    //void setSpinsChequer();
    void setSpinsSpecies(libconfig::Setting&);
    void setSpinsVampire(libconfig::Setting&);
    void setSpinsResText(libconfig::Setting&);
    void setSpinsGrid(libconfig::Setting&);
}
#endif /*_SPINS_H_*/
