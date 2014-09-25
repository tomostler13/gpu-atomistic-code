// File: anis.cpp
// Author: Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 25 Sep 2014 12:44:18
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/anis.h"
#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>
//#define FIXOUT(a,b) a.width(75);a << std::left << b;
//#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;
//#define SUCCESS(a) a << "Done" << std::endl;
#include "../inc/defines.h"

namespace anis
{
    Array<double> k1u;
    Array2D<double> k1udir;
}
