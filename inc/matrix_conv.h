// File: matrix_conv.h
// Author:Tom Ostler
// Created: 07 Oct 2014
// Last-modified: 08 Oct 2014 09:30:55
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _MATCONV_H_
#define _MATCONV_H_
namespace matconv
{
    void conv_intmat_to_dia(Array<int>&,Array4D<double>&,unsigned int&,unsigned int&,Array<double>&,Array<double>&,Array<double>&);
    void conv_intmat_to_dia(Array<int>&,Array<int>&,Array4D<double>&,unsigned int&,unsigned int&,unsigned int&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&,Array<double>&);
}
#endif /*_MATCONV_H_*/
