// File: spins.h
// Author:Tom Ostler
// Last-modified: 05 Dec 2012 12:06:03
#include "../inc/array4d.h"
#include "../inc/array3d.h"
#include "../inc/array.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include "../inc/config.h"
#ifndef _SPINS_H_
#define _SPINS_H_
namespace spins
{
	extern Array<double> sx;
	extern Array<double> sy;
	extern Array<double> sz;
	extern Array<double> esx;
	extern Array<double> esy;
	extern Array<double> esz;

    extern double mx;
    extern double my;
    extern double mz;
    extern double modm;

	extern Array3D<fftw_complex> csx;
	extern Array3D<fftw_complex> csy;
	extern Array3D<fftw_complex> csz;

	extern Array3D<fftw_complex> rsx;
	extern Array3D<fftw_complex> rsy;
	extern Array3D<fftw_complex> rsz;


	extern std::string sc;
	extern std::ifstream sfs;
	extern int sn[];
	extern bool si;
	extern fftw_plan SP;
	void initSpins(int argc,char *argv[]);
    void forwardTransformSpinArrays();
    void destroySpinPlans();
}
#endif /*_SPINS_H_*/
