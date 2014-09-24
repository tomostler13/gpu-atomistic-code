// File: config_glob.cpp
// Author:Tom Ostler
// Last-modified: 23 Sep 2014 14:44:51
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdarg>
#include "../inc/error.h"
#include "../inc/random.h"
#include "../inc/util.h"
#include "../inc/config.h"
#include <cassert>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
//File to initialise variables to tidy up main part of config.cpp
namespace config
{
    libconfig::Config cfg;
    bool lcf=false;
    bool inc_dip=false;
    unsigned int seed=0;
    std::ofstream Info,Log;
}
