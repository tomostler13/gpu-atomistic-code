// File: exch.h
// Author:Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 18 Jan 2013 14:10:19
#ifndef _EXCH_H_
#define _EXCH_H_
#include "../inc/arrays.h"
namespace exch
{
    extern Array<unsigned int> numint;
    extern Array2D<double> exchvec;
    extern Array2D<unsigned int> kvec;
    extern Array3D<double> J;
    extern unsigned int num_shells;
    void initExch(int argc,char *argv[]);
}
#endif /*_EXCH_H_*/
