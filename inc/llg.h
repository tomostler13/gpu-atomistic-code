// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 06 Jun 2015 13:56:19
#include "../inc/arrays.h"
#include <string>
#ifndef _LLG_H_
#define _LLG_H_
namespace llg
{
	extern double applied[],T,dt,rdt,gyro,muB;
    extern Array<double> llgpf;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&);
    extern std::string scm;
}
#endif /*_LLG_H_*/
