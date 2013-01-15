// File: config.h
// Author:Tom Ostler
// Last-modified: 05 Dec 2012 16:04:42
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
    extern std::ofstream Info;
    void initConfig(int argc,char *argv[]);
    void printline(std::ofstream&);
    std::string isTF(bool);
}
#endif /*_CONFIG_H_*/
