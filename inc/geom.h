// File: geom.h
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 17 Jan 2013 14:41:00
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _GEOM_H_
#define _GEOM_H_
namespace geom
{
    extern unsigned int dim[],nauc,zpdim[],nspins,zps;
    void initGeom(int argc,char *argv[]);
    extern Array2D<double> L,Linv;
    extern Array<double> abc;
    extern Array<unsigned int> Nk;


    extern unitCellMembers ucm;

}
#endif /*_GEOM_H_*/
