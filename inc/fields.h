// File: fields.h
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 15 May 2023 03:44:15 PM
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
    extern Array2D<double> sofield;
    extern Array5D<fftw_complex> Hk;
    extern Array5D<fftw_complex> Hr;
    extern Array4D<fftw_complex> dipHr,dipHk,hHk,hHr;
    extern Array<double> Hstagx,Hstagy,Hstagz,Hx,Hy,Hz,Hthx,Hthy,Hthz,HDemagx,HDemagy,HDemagz,H4sx,H4sy,H4sz;
    void initFields();
    void bfdip();
    void ftdip();
    void eftdip();
    void dipbfdip();
    void dipftdip();
    void dipeftdip();
    void setStagZero();
    void setStagField();
    void incStagFields(double);

}
#endif /*_FIELDS_H_*/
