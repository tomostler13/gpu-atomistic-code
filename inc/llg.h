// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 09 Apr 2013 21:23:01
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
    extern bool rscf;
    extern std::string osk;
    extern bool osHapp;
    extern bool osTemp;

    extern std::string rscfstr;
    extern std::ofstream rscfs;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&,Array<double>&);
    void integrate(unsigned int&);
    void integrate(unsigned int&,Array<double>&,Array<unsigned int>&,Array<unsigned int>&);
    void integrate(unsigned int&,Array<unsigned int>&,Array<unsigned int>&);
}
#endif /*_LLG_H_*/
