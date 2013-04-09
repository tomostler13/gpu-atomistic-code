// File: exch.h
// Author:Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 09 Apr 2013 15:11:53
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
    extern Array<unsigned int> xadj,adjncy;
    extern Array<double> Jxx,Jyy,Jzz;
    extern bool pbc [];
}
#endif /*_EXCH_H_*/
