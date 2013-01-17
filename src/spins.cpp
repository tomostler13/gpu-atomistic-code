// File: spins.cpp
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 17 Jan 2013 20:14:29
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/mat.h"
#include "../inc/intmat.h"
namespace spins
{
    Array3D<fftw_complex> Skx,Sky,Skz;
    Array<double> Sx,Sy,Sz;
    fftw_plan SxP,SyP,SzP;
    void initSpins(int argc,char *argv[])
    {
        Skx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Sky.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Skz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Sx.resize(geom::nspins);
        Sy.resize(geom::nspins);
        Sz.resize(geom::nspins);
        Sx.IFill(0);
        Sy.IFill(0);
        Sz.IFill(1);
        //forward transform of spin arrays
        SxP=fftw_plan_dft_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Skx.ptr(),Skx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        SyP=fftw_plan_dft_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Sky.ptr(),Sky.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        SzP=fftw_plan_dft_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Skz.ptr(),Skz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);

    }
    void FFTForward()
    {
        Skx.IFill(0);
        Sky.IFill(0);
        Skz.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Skx(lc[0],lc[1],lc[2])[0]=Sx[i];
            Sky(lc[0],lc[1],lc[2])[0]=Sy[i];
            Skz(lc[0],lc[1],lc[2])[0]=Sz[i];
        }
        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
    }
}
