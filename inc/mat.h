// File: mat.h
// Author:Tom Ostler
// Last-modified: 19 Dec 2012 10:12:10
#include <string>
#include "../inc/array.h"
#include "../inc/array2d.h"
#ifndef _MAT_H_
#define _MAT_H_
namespace mat
{
    extern bool mi;
	extern double Ms;
	extern double lambda;
	extern double Tc;
	extern double gamma;
    extern double K;
    extern double chippf;
    extern std::string materialid;
    extern unsigned int nlat;
    extern Array<double> mu;
    extern Array<double> conc;
	void initMat(int argc,char *argv[]);
}
#endif /*_MAT_H_*/
