// File: mat.h
// Author:Tom Ostler
// Last-modified: 22 Apr 2013 15:33:03
#include <string>
#include "../inc/array.h"
#include "../inc/array2d.h"
#ifndef _MAT_H_
#define _MAT_H_
namespace mat
{
    extern unsigned int nspec;
    extern double gamma,muB,ilambda;
    extern Array<double> sigma,mu,mustore,lambda;
    extern std::string place;
    extern Array<unsigned int> speclist;
    void initMat(int argc,char *argv[]);
}
#endif /*_MAT_H_*/
