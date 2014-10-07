// File: matrix_conv.h
// Author:Tom Ostler
// Created: 07 Oct 2014
// Last-modified: 07 Oct 2014 11:06:24
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _MATCONV_H_
#define _MATCONV_H_
namespace matconv
{
    void dia_offsets(Array<int>&,Array4D<double>&,unsigned int&,unsigned int&);
}
#endif /*_MATCONV_H_*/
