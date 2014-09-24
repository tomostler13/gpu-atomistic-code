// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 24 Sep 2014 11:34:23
#include "../inc/arrays.h"
#ifndef _LLG_H_
#define _LLG_H_
namespace llg
{
	extern double applied[],T,dt,rdt,gyro,muB;
    extern Array<double> llgpf;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&);
}
#endif /*_LLG_H_*/
