// File: config.h
// Author:Tom Ostler
// Last-modified: 09 Apr 2013 15:43:45
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include "../inc/error.h"
#ifndef _CONFIG_H_
#define _CONFIG_H_
namespace config
{
    extern libconfig::Config cfg;
    extern unsigned int seed;
    extern bool lcf;
    extern bool incdip;
    extern std::ofstream Info;
    void initConfig(int argc,char *argv[]);
    void printline(std::ofstream&);
    std::string isTF(bool);
    extern bool useintmat;
}
#endif /*_CONFIG_H_*/
