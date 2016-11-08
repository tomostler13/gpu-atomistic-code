// File: rscf.h
// Author:Tom Ostler
// Created: 23 Oct 2015
// Last-modified: 18 Oct 2016 13:26:44
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#include <fstream>
#ifndef _RSCF_H_
#define _RSCF_H_
namespace rscf
{
    extern bool ccf;
    void initRSCF();
    void calcRSCF(unsigned int);
    std::string conv(unsigned int);
    void outputRSCF(std::ofstream&,Array3D<double>&,unsigned int);
    void cleanup();
}
#endif /*_RSCF_H_*/
