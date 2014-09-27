// File: fields.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 26 Sep 2014 16:14:10
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../inc/arrays.h"
#ifndef _FIELDS_H_
#define _FIELDS_H_
namespace fields
{
    extern Array5D<fftw_complex> Hk;
    extern Array5D<fftw_complex> Hr;
    extern Array<double> Hx,Hy,Hz,Hthx,Hthy,Hthz;
    void initFields(int argc,char *argv[]);
    void bfdip();
    void ftdip();
    void eftdip();
}
#endif /*_FIELDS_H_*/
