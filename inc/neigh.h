// File: neigh.h
// Author:Tom Ostler
// Last-modified: 03 Dec 2012 12:13:48
#include <cstdlib>
#include "../inc/array.h"
#ifndef _NEIGH_H_
#define _NEIGH_H_
namespace neigh
{
    extern bool pbc[];
    extern bool ni;
    extern Array<unsigned int> xadj;
    extern Array<unsigned int> adjncy;
    extern Array<double> surfArea;

    void initNeigh(int argc,char *argv[]);
    int bc(int,int);
}
#endif /*_NEIGH_H_*/
