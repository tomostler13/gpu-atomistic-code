// File: llg.h
// Author:Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 21 Jan 2013 15:17:31
#ifndef _LLG_H_
#define _LLG_H_
namespace llg
{
    extern double dt,rdt,applied[],T,llgpf;
    void initLLG(int argc,char *argv[]);
    void llgCPU(unsigned int);
}
#endif /*_LLG_H_*/
