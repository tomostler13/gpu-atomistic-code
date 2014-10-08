// File: matrix_mul_cpu.h
// Author:Tom Ostler
// Created: 08 Oct 2014
// Last-modified: 08 Oct 2014 16:20:25
#include "../inc/arrays.h"
#include <string>
#ifndef _MATMUL_H_
#define _MATMUL_H_
namespace matmul
{
    void spmv_dia_diag(unsigned int,int,Array<int>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&);
    void spmv_dia_offdiag(unsigned int,int,Array<int>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&);
}
#endif /*_MATMUL_H_*/
