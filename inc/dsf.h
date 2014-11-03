// File: dsf.h
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 03 Nov 2014 13:17:14
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _DSF_H_
#define _DSF_H_
namespace dsf
{
    extern Array2D<double> Ry,Rz,uo;
    extern Array2D<unsigned int> symdir;
    extern Array<double> qa;
    extern bool cdsf;
    extern double et,rt,T;
    extern unsigned int dsfupdate,nsd,ets,rts;
    void initDSF(int argc,char *argv[]);
    void readDSFParam(int argc,char *argv[]);
}
#endif /*_DSF_H_*/
