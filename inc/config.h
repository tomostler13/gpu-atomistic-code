// File: config.h
// Author:Tom Ostler
// Last-modified: 03 Oct 2014 10:10:53
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
    //Config initialised
    extern bool lcf;
    //include dipolar fields?
    extern bool inc_dip;
    extern std::string intmeth;
    extern unsigned int intm;
    extern std::ofstream Info,Log;
    void initConfig(int argc,char *argv[]);
    void printline(std::ofstream&);
    void openLogFile();
    std::string isTF(bool);
}
#endif /*_CONFIG_H_*/
