// File: exch.h
// Author:Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 05 Jun 2015 17:24:00
#ifndef _EXCH_H_
#define _EXCH_H_
#include <string>
#include "../inc/arrays.h"
namespace exch
{
    extern Array<int> diagoffset,offdiagoffset,checkdiag;
    extern Array<unsigned int> xadj,adjncy,offdiagxadj,offdiagadjncy,xadj_j,xadj_k,xadj_l,adjncy_j,adjncy_k,adjncy_l;
    extern Array2D<unsigned int> numquart;
    extern Array<double> dataxx,dataxy,dataxz,datayx,datayy,datayz,datazx,datazy,datazz;
    extern Array<double> evs;
    extern Array2D<unsigned int> shell_list;
    extern Array3D<unsigned int> numint;
    extern Array4D<double> exchvec;
    extern Array5D<int> fsqi;
    extern Array5D<double> J,fsq;
    extern unsigned int num_shells,diagnumdiag,offdiagnumdiag,max_shells,max_int,max_4s;
    extern bool outputJ,oem,rem,cutexch,inc4spin;
    extern double rcut;
    extern std::string readMethod,readFile,method,enerType;
    void initExch(int argc,char *argv[]);
    void read4spin();
    void readGlobalExch(int argc,char *argv[]);

    extern libconfig::Config exchcfg;
    void permute(int argc,char *argv[]);
    void direct(int argc,char *argv[]);
    void mapint(int argc,char *argv[]);
    void get_exch_permute(int argc,char *argv[]);
    void get_exch_direct(int argc,char *argv[]);
    void get_exch_mapint(int argc,char *argv[]);
    void hybrid();
    void fft();
    void directfft();
    void directhybrid();
    void exchpermute();
    void exchdirect();
    void exchmapint();
    void setup4SpinCSR();
}
#endif /*_EXCH_H_*/
