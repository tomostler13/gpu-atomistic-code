// File: geom.h
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 08 Jan 2018 20:47:10
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <string>
#ifndef _GEOM_H_
#define _GEOM_H_
namespace geom
{
    extern unsigned int dim[],zpdim[],nspins,zps,cplxdim,czps,nms,nauc;
    extern std::string place;
    void initGeom();
    void readconfig();
    extern bool Nkset,rprimset;
    extern Array2D<int> lu,zplu;
    extern Array2D<double> rprim;
    extern Array4D<int> coords;
    extern Array<double> abc,mu,gamma,lambda,llgpf,rx,ry,rz,sigma;
    extern Array<unsigned int> Nk,sublattice;
    extern bool logunit;
    extern Array4D<unsigned int> atnolu;

    extern unitCellMembers ucm;
}
#endif /*_GEOM_H_*/
