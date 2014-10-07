// File: config.h
// Author:Tom Ostler
// Last-modified: 07 Oct 2014 10:40:07
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
    extern bool lcf,offdiag;
    //include dipolar fields?
    extern bool inc_dip,pbc[];
    extern std::string exchmeth,dipmeth;
    extern unsigned int exchm,dipm;
    extern std::ofstream Info,Log;
    void initConfig(int argc,char *argv[]);
    void printline(std::ofstream&);
    void openLogFile();
    std::string isTF(bool);
}
#endif /*_CONFIG_H_*/
