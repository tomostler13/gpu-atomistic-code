// File: geom_glob.cpp
// Author:Tom Ostler
// Created: 26 July 2014
// Last-modified: 26 Jun 2014 10:11:28
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

//To keep the code tidy the declaration of the geom namespace variables are declared here
namespace geom
{
    //the number of unit cells
    unsigned int dim[3]={0,0,0};
    //zero pad size (in unit cells)
    unsigned int zpdim[3]={0,0,0};
    //For r2c transform the final dimension must be zpdim[2]/2+1
    unsigned int cplxdim=0;
    //the number of atoms in the unit cell
    unsigned int nauc=0;
    //maximum system size, number of spins, zeropad size
    //number of magnetic species
    unsigned int maxss=0,nspins=0,zps=0,czps=0,nms=0;
    //The unit vectors describing the lattice unit cell (unit vectors)
    Array2D<double> L,Linv;
    //The a,b and c values (i.e. lattice constants)
    Array<double> abc;
    //Number of K points
    Array<unsigned int> Nk;
    //lookup array. Give atom number and return coordinates
    Array2D<int> lu;
    //Coords array. Give coords and returns atom number. This coords array
    //is the size of the zero padded arrays and has places where atoms do
    //not exist
    std::string place;
    Array4D<int> coords;
    //instance of the unit cell
    unitCellMembers ucm;
}
