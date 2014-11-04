// File: dsf.h
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 04 Nov 2014 16:54:10
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _DSF_H_
#define _DSF_H_
namespace dsf
{
    extern Array2D<double> Ry,Rz,uo;
    extern Array2D<unsigned int> kpoints;
    extern Array<double> qa;
    extern bool cdsf;
    extern double et,rt,T;
    extern unsigned int dsfupdate,nk,ets,rts;
    void initDSF(int argc,char *argv[]);
    void readDSFParam(int argc,char *argv[]);
}
#endif /*_DSF_H_*/
