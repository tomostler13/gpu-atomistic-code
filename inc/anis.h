// File: anis.h
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 18 Mar 2014 17:01:25
#ifndef _ANIS_H_
#define _ANIS_H_
#include "../inc/arrays.h"
namespace anis
{
    void initAnis(int argc,char *argv[]);
    extern unsigned int nfou;
    extern Array<double> FirstOrderUniaxK;
    extern Array2D<double> FirstOrderUniaxDir;
}
#endif /*_ANIS_H_*/
