// File: rscf.h
// Author:Tom Ostler
// Created: 23 Oct 2015
// Last-modified: 21 Feb 2017 10:29:22
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#include <fstream>
#ifndef _RSCF_H_
#define _RSCF_H_
namespace rscf
{
    extern bool ccf;
    extern unsigned int upd;
    void initRSCF();
    void calcRSCF(unsigned int);
    std::string conv(unsigned int);
    void outputRSCF(std::ofstream&,Array3D<double>&,unsigned int);
    void outputRSCFNextIndex();
    void cleanup();
}
#endif /*_RSCF_H_*/
