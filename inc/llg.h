// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 03 Aug 2015 15:38:27
#include "../inc/arrays.h"
#include <string>
#ifndef _LLG_H_
#define _LLG_H_
namespace llg
{
	extern double applied[],T,dt,rdt,gyro,muB,kB;
    extern Array<double> llgpf,Ts,cps,dps;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&);
    extern std::string scm;
}
#endif /*_LLG_H_*/
