// File: LLBCPU.h
// Author:Tom Ostler
// Last-modified: 04 Dec 2012 17:21:09
#include "../inc/array2d.h"
#ifndef _LLBCPU_H_
#define _LLBCPU_H_
namespace llb
{
    extern bool LLBi;
    extern Array2D<double> fn;
    void initLLBCPU(int argc, char *argv[]);
    void LLBCPU(unsigned int&);
}
#endif /*_LLBCPU_H_*/
