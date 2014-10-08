// File: fields.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 08 Oct 2014 14:24:13
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
    extern Array4D<fftw_complex> dipHr,dipHk;
    extern Array<double> Hx,Hy,Hz,Hthx,Hthy,Hthz;
    void initFields(int argc,char *argv[]);
    void bfdip();
    void ftdip();
    void eftdip();
    void dipbfdip();
    void dipftdip();
    void dipeftdip();

}
#endif /*_FIELDS_H_*/
