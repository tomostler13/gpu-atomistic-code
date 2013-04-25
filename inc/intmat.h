// File: intmat.h
// Author:Tom Ostler
// Created: 16 Jan 2012
// Last-modified: 24 Apr 2013 22:10:34
#include "../inc/arrays.h"
#include <fftw3.h>
#ifndef _INTMAT_H_
#define _INTMAT_H_
namespace intmat
{
    extern Array3D<fftw_complex> Nxx,Nxy,Nxz,Nyx,Nyy,Nyz,Nzx,Nzy,Nzz;
    extern Array3D<double> Nrxx,Nrxy,Nrxz,Nryx,Nryy,Nryz,Nrzx,Nrzy,Nrzz;
    extern Array<unsigned int> zpsn;
    void initIntmat(int argc,char *argv[]);
    void fillIntmat();
    void fftIntmat();
}
#endif /*_INTMAT_H_*/
