// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 25 Mar 2013 17:31:02
#ifndef _LLG_H_
#define _LLG_H_
#include <string>
#include <sstream>
#include <fstream>
namespace llg
{
	extern double applied[],T,dt,rdt,llgpf;
    extern bool rscf;
    extern std::string rscfstr;
    extern std::ofstream rscfs;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&);
}
#endif /*_LLG_H_*/
