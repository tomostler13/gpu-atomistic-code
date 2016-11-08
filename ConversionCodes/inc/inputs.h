// File: inputs.h
// Author:Tom Ostler
// Created: 25 Jan 2016
// Last-modified: 25 Feb 2016 22:20:02
#include "../../inc/arrays.h"
#include <string>
#include <fstream>
#ifndef _INPUTS_H_
#define _INPUTS_H_
namespace inputs
{
    extern std::ofstream Info;
    extern unsigned int nspins,lt,ut,dt;
    extern unsigned int dimin[],dimout[],Nk[];
    extern double abc[];
    extern std::string fb;
    void readcff(int argc,char *argv[]);
}
#endif /*_INPUTS_H_*/
