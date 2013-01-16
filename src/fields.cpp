// File: fields.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 16 Jan 2013 17:59:22
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/mat.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace fields
{
    Array3D<fftw_complex> Hkx,Hky,Hkz;
    Array<double> Hx,Hy,Hz;
}
