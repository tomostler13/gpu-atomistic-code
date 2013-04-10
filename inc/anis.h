// File: anis.h
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 10 Apr 2013 13:47:00
#ifndef _ANIS_H_
#define _ANIS_H_
#include "../inc/arrays.h"
namespace anis
{
    void initAnis(int argc,char *argv[]);
    extern Array2D<double> dT;
    extern double uniaxial_unit[];
}
#endif /*_ANIS_H_*/
