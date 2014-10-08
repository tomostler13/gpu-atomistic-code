// File: intmat.h
// Author:Tom Ostler
// Created: 16 Jan 2012
// Last-modified: 08 Oct 2014 13:54:51
#include "../inc/arrays.h"
#ifndef _INTMAT_H_
#define _INTMAT_H_
namespace intmat
{
    extern Array7D<fftw_complex> Nkab;
    extern Array7D<fftw_complex> Nrab;
    extern Array5D<fftw_complex> dipNrab,dipNkab;
    extern Array<unsigned int> zpsn;
    void initIntmat(int argc,char *argv[]);
    void initDipIntmat(int argc,char *argv[]);
    void fillIntmat();
    void fftIntmat();
    void fftDipIntmat();
    void fillDipIntmat();
}
#endif /*_INTMAT_H_*/
