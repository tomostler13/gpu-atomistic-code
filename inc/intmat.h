// File: intmat.h
// Author:Tom Ostler
// Created: 16 Jan 2012
// Last-modified: 24 Sep 2014 14:35:10
#include "../inc/arrays.h"
#ifndef _INTMAT_H_
#define _INTMAT_H_
namespace intmat
{
    extern Array7D<fftw_complex> Nkab;
    extern Array7D<double> Nrab;
    extern Array<unsigned int> zpsn;
    void initIntmat(int argc,char *argv[]);
    void fillIntmat();
    void fftIntmat();
}
#endif /*_INTMAT_H_*/
