// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 18 May 2023 10:14:50 AM
#include "../inc/arrays.h"
#include <string>
#ifndef _LLG_H_
#define _LLG_H_
namespace llg
{
	extern double applied[],T,dt,rdt,gyro,muB,kB;
    extern Array2D<double> optrans;
    extern unsigned int trans,intscheme;
    extern bool exstr;
    extern Array<double> llgpf,Ts,cps,dps;
	void initLLG();
    void integrate(unsigned int&);
    extern std::string scm,intschemestr;
}
#endif /*_LLG_H_*/
