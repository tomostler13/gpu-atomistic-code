// File: exch.h
// Author:Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 05 Dec 2014 11:18:36
#ifndef _EXCH_H_
#define _EXCH_H_
#include <string>
#include "../inc/arrays.h"
namespace exch
{
    extern Array<int> diagoffset,offdiagoffset,checkdiag;
    extern Array<unsigned int> xadj,adjncy,offdiagxadj,offdiagadjncy;
    extern Array<double> dataxx,dataxy,dataxz,datayx,datayy,datayz,datazx,datazy,datazz;
    extern Array<double> evs;
    extern Array2D<unsigned int> shell_list;
    extern Array3D<unsigned int> numint;
    extern Array4D<double> exchvec;
    extern Array5D<double> J;
    extern unsigned int num_shells,diagnumdiag,offdiagnumdiag,max_shells;
    extern bool outputJ,oem,rem;
    extern std::string readMethod,readFile,method,enerType;
    void initExch(int argc,char *argv[]);
    void readGlobalExch(int argc,char *argv[]);

    extern libconfig::Config exchcfg;
    void permute();
    void direct();
}
#endif /*_EXCH_H_*/
