// File: mat.h
// Author:Tom Ostler
// Last-modified: 23 Sep 2014 10:56:04
#include <string>
#include "../inc/array.h"
#include "../inc/array2d.h"
#ifndef _MAT_H_
#define _MAT_H_
namespace mat
{
    extern Array<double> lambda,gamma,mu,sigma;
    extern double muB,gyro;
    void initMat(int argc,char *argv[]);
}
#endif /*_MAT_H_*/
