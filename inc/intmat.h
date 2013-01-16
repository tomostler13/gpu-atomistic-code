// File: intmat.h
// Author:Tom Ostler
// Created: 16 Jan 2012
// Last-modified: 16 Jan 2013 18:02:48
#include "../inc/array3d.h"
#ifndef _INTMAT_H_
#define _INTMAT_H_
namespace intmat
{
    extern Array3D<fftw_complex> Nxx,Nxy,Nxz,Nyx,Nyy,Nyz,Nzx,Nzy,Nzz;
    extern Array<unsigned int> zpsn;
    void initIntmat(int argc,char *argv[]);
    void fillIntmat();
}
#endif /*_INTMAT_H_*/
