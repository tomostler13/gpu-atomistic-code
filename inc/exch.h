// File: exch.h
// Author:Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 10 Oct 2014 14:56:35
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
    extern Array4D<unsigned int> kvec;
    extern Array5D<double> J;
    extern unsigned int num_shells,diagnumdiag,offdiagnumdiag;
    void initExch(int argc,char *argv[]);
}
#endif /*_EXCH_H_*/
