// File: llg.cpp
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 22 Jan 2013 15:38:01
#include "../inc/LLG.h"
#include "../inc/LLGCPU.h"
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace LLG
{
    void integrate(unsigned int& t)
    {
        #ifdef CUDA
        cullg::llgGPU(t);
        #else
        llbg::llgCPU(t);
        #endif
    }
}
