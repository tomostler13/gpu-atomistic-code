// File: matrix_conv.h
// Author:Tom Ostler
// Created: 07 Oct 2014
// Last-modified: 07 Oct 2014 13:14:42
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include "../inc/exch.h"
#include <string>
#ifndef _MATCONV_H_
#define _MATCONV_H_
namespace matconv
{
    void conv_intmat_to_dia(Array<int>&,Array4D<double>&,unsigned int&,unsigned int&,Array<double>&);
}
#endif /*_MATCONV_H_*/
