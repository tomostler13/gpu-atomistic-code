// File: fields.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 25 Mar 2013 12:18:24
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
    extern Array3D<double> Hrx,Hry,Hrz;
    extern Array<double> Hx,Hy,Hz,Hthx,Hthy,Hthz;
    void initFields(int argc,char *argv[]);
    void bfdip();
    void ftdip();
    void eftdip();
}
#endif /*_FIELDS_H_*/
