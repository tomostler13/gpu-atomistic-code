// File: geom.h
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 24 Sep 2014 09:22:40
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _GEOM_H_
#define _GEOM_H_
namespace geom
{
    extern unsigned int dim[],zpdim[],nspins,zps,cplxdim,czps,nms;
    extern std::string place;
    void initGeom(int argc,char *argv[]);
    void readconfig(int argc,char *argv[]);
    extern Array2D<int> lu,zplu;
    extern Array4D<int> coords;
    extern Array<double> abc,gamma,lambda,llgpf,rx,ry,rz,sublattice;
    extern Array<unsigned int> Nk;

    extern unitCellMembers ucm;
}
#endif /*_GEOM_H_*/
