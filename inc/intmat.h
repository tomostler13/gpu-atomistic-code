// File: intmat.h
// Author:Tom Ostler
// Created: 16 Jan 2012
// Last-modified: 18 Oct 2016 12:52:43
#include "../inc/arrays.h"
#ifndef _INTMAT_H_
#define _INTMAT_H_
namespace intmat
{
    extern Array7D<fftw_complex> Nkab;
    extern Array7D<fftw_complex> Nrab;
    extern Array5D<fftw_complex> dipNrab,dipNkab,hNrab,hNkab;
    extern Array<unsigned int> zpsn;
    void initIntmat();
    void initDipIntmat();
    void fillIntmat();
    void fftIntmat();
    void fftDipIntmat();
    void fillDipIntmat();
}
#endif /*_INTMAT_H_*/
