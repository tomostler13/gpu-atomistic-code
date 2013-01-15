// File: LLB.cpp
// Author:Tom Ostler
// Last-modified: 04 Jan 2013 14:36:23
#include "../inc/LLB.h"
#include "../inc/LLBCPU.h"
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace LLB
{
    void integrate(unsigned int& t)
    {
        #ifdef CUDA
        assert(cullb::cullbinit);
        cullb::LLBGPU(t);
        #else
        llb::LLBCPU(t);
        #endif
    }
}
