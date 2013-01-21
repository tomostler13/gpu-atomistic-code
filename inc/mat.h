// File: mat.h
// Author:Tom Ostler
// Last-modified: 21 Jan 2013 15:19:17
#include <string>
#include "../inc/array.h"
#include "../inc/array2d.h"
#ifndef _MAT_H_
#define _MAT_H_
namespace mat
{
    extern double lambda,gamma,muB,mu,sigma,sigma;
    void initMat(int argc,char *argv[]);
}
#endif /*_MAT_H_*/
