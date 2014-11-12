// File: exch.h
// Author:Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 11 Nov 2014 13:58:00
#ifndef _EXCH_H_
#define _EXCH_H_
#include "../inc/arrays.h"
namespace exch
{
    extern Array<int> diagoffset,offdiagoffset;
    extern Array<unsigned int> xadj,adjncy,offdiagxadj,offdiagadjncy;
    extern Array<double> dataxx,dataxy,dataxz,datayx,datayy,datayz,datazx,datazy,datazz;
    extern Array3D<unsigned int> numint;
    extern Array4D<double> exchvec;
    extern Array5D<double> J;
    extern unsigned int num_shells,diagnumdiag,offdiagnumdiag;
    extern bool outputJ;
    void initExch(int argc,char *argv[]);
}
#endif /*_EXCH_H_*/
