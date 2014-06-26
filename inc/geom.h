// File: geom.h
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 26 Jun 2014 10:13:29
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _GEOM_H_
#define _GEOM_H_
namespace geom
{
    extern unsigned int dim[],nauc,zpdim[],nspins,zps,cplxdim,czps,nms,maxss;
    extern std::string place;
    void initGeom(int argc,char *argv[]);
    extern Array2D<double> L,Linv;
    extern Array2D<int> lu,zplu;
    extern Array4D<int> coords;
    extern Array<double> abc;
    extern Array<unsigned int> Nk;

    extern unitCellMembers ucm;
}
#endif /*_GEOM_H_*/
