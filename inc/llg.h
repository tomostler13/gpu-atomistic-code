// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 17 Apr 2013 11:44:29
#ifndef _LLG_H_
#define _LLG_H_
#include "../inc/array.h"
#include <string>
#include <sstream>
#include <fstream>
namespace llg
{
	extern double applied[],T,dt,rdt,llgpf;
    extern Array<double> osT;
    extern bool rscf,ssf;
    extern std::string osk;
    extern bool osHapp;
    extern bool osTemp;

    extern std::string rscfstr,ssffs;
    extern std::ofstream rscfs,ssffstr;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&);
}
#endif /*_LLG_H_*/
