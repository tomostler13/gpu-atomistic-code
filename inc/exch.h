// File: exch.h
// Author:Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 07 Oct 2014 09:27:11
#ifndef _EXCH_H_
#define _EXCH_H_
#include "../inc/arrays.h"
namespace exch
{
    extern Array3D<unsigned int> numint;
    extern Array4D<double> exchvec;
    extern Array4D<unsigned int> kvec;
    extern Array5D<double> J;
    extern Array4D<double> JMat;
    extern unsigned int num_shells;
    void initExch(int argc,char *argv[]);
}
#endif /*_EXCH_H_*/
