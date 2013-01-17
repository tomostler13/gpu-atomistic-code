// File: fields.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 17 Jan 2013 14:47:18
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
    extern Array3D<fftw_complex> Hkx,Hky,Hkz;
    extern Array<double> Hx,Hy,Hz;
    void initFields(int argc,char *argv[]);
    void bfdip();
    void ftdip();
}
#endif /*_FIELDS_H_*/
