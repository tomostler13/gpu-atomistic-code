// File: anis.h
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 12 Apr 2023 04:47:01 PM
#ifndef _ANIS_H_
#define _ANIS_H_
#include "../inc/arrays.h"
namespace anis
{
    extern Array<double> k1u;
    extern Array2D<double> k1udir;
    extern double k2par,k2perp,k4par,k4perp;
    extern Array<double> k2perpdir,k2pardir;
    extern Array2D<double> k4perpdirs,k4pardirs;
    void readGlobalAnis();
}
#endif /*_ANIS_H_*/
