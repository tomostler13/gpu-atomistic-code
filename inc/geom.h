// File: geom.h
// Author:Tom Ostler
// Last-modified: 27 Dec 2012 14:34:13
#include "../inc/array3d.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include <string>
#ifndef _GEOM_H_
#define _GEOM_H_
namespace geom
{
    extern double gs[];
    extern unsigned int dim[];
    extern int zpdim[];
	extern unsigned int cplxzpdim;
    extern int zps;
    extern unsigned int ss;
    extern unsigned int maxss;
	extern bool gi;
    extern std::string structure;
    extern double gsV;
    extern Array2D<unsigned int> lu;
    extern Array3D<int> coords;

    void initGeom(int argc,char *argv[]);
}
#endif /*_GEOM_H_*/
