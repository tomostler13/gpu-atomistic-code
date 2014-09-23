// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 23 Sep 2014 10:44:31
#include "../inc/arrays.h"
#ifndef _LLG_H_
#define _LLG_H_
namespace llg
{
	extern double applied[],T,dt,rdt;
    extern Array<double> llgpf;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&);
}
#endif /*_LLG_H_*/
