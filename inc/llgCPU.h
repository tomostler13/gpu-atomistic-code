// File: llg.h
// Author:Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 09 Apr 2013 21:29:27
#include "../inc/array.h"
#ifndef _LLGCPU_H_
#define _LLGCPU_H_
namespace llgCPU
{
    void initLLG(int argc,char *argv[]);
    void llgCPU(unsigned int,Array<double>&);
    void llgCPU(unsigned int);
    void llgCPU(unsigned int,Array<double>&,Array<unsigned int>&,Array<unsigned int>&);
    void llgCPU(unsigned int,Array<unsigned int>&,Array<unsigned int>&);

}
#endif /*_LLG_H_*/
