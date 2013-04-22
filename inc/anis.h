// File: anis.h
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 22 Apr 2013 14:17:42
#ifndef _ANIS_H_
#define _ANIS_H_
#include "../inc/arrays.h"
namespace anis
{
    void initAnis(int argc,char *argv[]);
    extern Array3D<double> dT;
    extern Array2D<double> uniaxial_unit;
}
#endif /*_ANIS_H_*/
