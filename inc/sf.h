// File: sf.h
// Note: originally dsf.h
// Author:Tom Ostler
// Created: 2 Nov 2014
// Last-modified: 18 Oct 2016 13:26:53
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
    void initSF();
    void readSFParam();
}
#endif /*_SF_H_*/
