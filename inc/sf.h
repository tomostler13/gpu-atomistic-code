// File: sf.h
// Note: originally dsf.h
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 06 Nov 2014 12:01:44
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _SF_H_
#define _SF_H_
namespace sf
{
    extern Array2D<double> Ry,Rz,uo;
    extern Array2D<unsigned int> kpoints;
    extern Array<double> qa;
    extern bool csf;
    extern unsigned int sfupdate,nk,ets,rts;
    void initSF(int argc,char *argv[]);
    void readSFParam(int argc,char *argv[]);
}
#endif /*_SF_H_*/
