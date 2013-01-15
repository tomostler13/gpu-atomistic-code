// File: intmat.h
// Author:Tom Ostler
// Last-modified: 27 Dec 2012 11:53:24
#include "../inc/array3d.h"
#ifndef _INTMAT_H_
#define _INTMAT_H_
namespace intmat
{
	extern Array3D<fftw_complex> Nxx;
 	extern Array3D<fftw_complex> Nxy;
  	extern Array3D<fftw_complex> Nxz;
 	extern Array3D<fftw_complex> Nyx;
	extern Array3D<fftw_complex> Nyy;
	extern Array3D<fftw_complex> Nyz;
	extern Array3D<fftw_complex> Nzx;
	extern Array3D<fftw_complex> Nzy;
	extern Array3D<fftw_complex> Nzz;
    extern bool imi;
    extern bool imft;

	void initIntmat(int argc,char *argv[]);
	void fillIntmat();
	void fftIntmat();
}
#endif /*_INTMAT_H_*/
