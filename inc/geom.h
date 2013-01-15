// File: geom.h
// Author:Tom Ostler
// Last-modified: 15 Jan 2013 17:39:31
#include "../inc/array3d.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include <string>
#ifndef _GEOM_H_
#define _GEOM_H_
namespace geom
{
    extern unsigned int dim[];
    extern int zpdim[];
    extern int zps;
    extern unsigned int ss;
    extern unsigned int maxss;
	extern bool gi;
    extern std::string structure;
    extern unsigned int maxss;
    extern unsigned int nauc;
    extern Array2D<double> atompos;
    extern Array2D<unsigned int> lu;
    extern Array3D<int> coords;
    extern Array<unsigned int> scalecoords;

    void initGeom(int argc,char *argv[]);
}
#endif /*_GEOM_H_*/
