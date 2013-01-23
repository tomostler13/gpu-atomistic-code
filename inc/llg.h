// File: llg.h
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 22 Jan 2013 16:17:35
#ifndef _LLG_H_
#define _LLG_H_
namespace llg
{
	extern double applied[],T,dt,rdt,llgpf;
	void initLLG(int argc,char *argv[]);
    void integrate(unsigned int&);
}
#endif /*_LLG_H_*/
