// File: tdp.h
// Author:Tom Ostler
// Last-modified: 31 Dec 2012 15:51:11
#include <string>
#include <cstdlib>
#include "../inc/array.h"
#ifndef _TDP_H_
#define _TDP_H_
namespace tdp
{
    extern bool tdpi;
    extern double uniformtemp;
    extern std::string toh;
    extern Array<double> systemp;
    extern Array<double> syschipar;
    extern Array<double> syschiperp;
    extern Array<double> sysexchstiff;
    extern Array<double> sysme;
    extern Array<double> sysalphaperp;
    extern Array<double> sysalphapar;
    extern Array<double> sysW1pf;
    extern Array<double> sysW2pf;
    extern Array<double> sysSigma;
    extern std::string chiperptype;
    extern std::string chipartype;
    extern std::string exchstifftype;
    extern std::string metype;
    extern unsigned int chiperpfunc;
    extern unsigned int chiparfunc;
    extern unsigned int mefunc;
    extern unsigned int exchstifffunc;
    extern bool tss;
    extern double tssf;
    extern void (*Tcalc)();
    extern void (*chiperpfp)();
    extern void (*chiparfp)();
    extern void (*exchstifffp)();
    extern void (*mefp)();
    extern double chiperpsingle(double);
    void inittdp(int argc,char *argv[]);
    void Tcalc0u();
    void chiperp0();
    void chiperpmf();
    void chipar0();
    void chiparmf();
    void exchstiff0();
    void exchstiffNULL();
    void me0();
    void memf();
    void calcalphas();
}
#endif /*_TDP_H_*/
