// File: fields.h
// Author:Tom Ostler
// Last-modified: 10 Dec 2012 20:21:40
#include <fftw3.h>
#include "../inc/array3d.h"
#include "../inc/array.h"
#include <cstdlib>
#include <string>
#ifndef _FIELDS_H_
#define _FIELDS_H_
namespace fields
{
	extern Array3D<fftw_complex> Hkx;
	extern Array3D<fftw_complex> Hky;
	extern Array3D<fftw_complex> Hkz;

    extern Array<double> Hx;
    extern Array<double> Hy;
    extern Array<double> Hz;

    extern Array<double> GW1x;
    extern Array<double> GW1y;
    extern Array<double> GW1z;

    extern Array<double> GW2x;
    extern Array<double> GW2y;
    extern Array<double> GW2z;

    extern bool checkthermal;
    extern bool checkdipolar;
    extern bool checkexchange;

    extern unsigned int dfu;
    extern unsigned int anist;
    extern void (*anisfp)(double *,const double *,unsigned int);

    extern unsigned int ff;
    extern void (*appliedfp)(double *);
    extern double fH[];

	extern bool fi;
	extern std::string dfc;
	extern fftw_plan HxP,HyP,HzP;

	void initFields(int argc,char *argv[]);
    void transformFieldsBack();
    void destroyFieldPlans();
    void bfdip();
    void anis0(double *,const double *,unsigned int);
    void anis1u(double *,const double *,unsigned int);
    void applied0(double *);
}
#endif /*_FIELDS_H_*/
