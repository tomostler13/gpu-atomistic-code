// File: maths.cpp
// Author:Tom Ostler
// Last-modified: 05 Dec 2012 16:07:06
#include "../inc/maths.h"
#include "../inc/random.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/geom.h"
#include "../inc/intmat.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <fftw3.h>
namespace maths
{
	void sphereRandom(double& x,double& y,double& z)
	{
		double v1=0,v2=0,s=2.0,ss=0.0;
        unsigned int counter=0;
		while(s>1.0)
		{
			v1=2.0*Random::rand()-1.0;
			v2=2.0*Random::rand()-1.0;
			s=v1*v1+v2*v2;
            counter++;
		}
		ss=sqrt(1.0-s);
		x=2.0*v1*ss;
		y=2.0*v2*ss;
		z=1.0-2.0*s;
	}

	void cpuConvFourier()
	{
		assert(fields::fi);
		assert(spins::si);
		assert(geom::gi);
        spins::forwardTransformSpinArrays();
        fields::Hkx.IFill(0);
        fields::Hky.IFill(0);
        fields::Hkz.IFill(0);

		//perform convolution in fourier space
		for(unsigned int i = 0 ; i < geom::zpdim[0] ; i++)
		{
			for(unsigned int j = 0 ; j < geom::zpdim[1] ; j++)
			{
				for(unsigned int k = 0 ; k < geom::zpdim[2] ; k++)
				{

                    fields::Hkx(i,j,k)[0]=intmat::Nxx(i,j,k)[0]*spins::csx(i,j,k)[0]-intmat::Nxx(i,j,k)[1]*spins::csx(i,j,k)[1]
                                         +intmat::Nxy(i,j,k)[0]*spins::csy(i,j,k)[0]-intmat::Nxy(i,j,k)[1]*spins::csy(i,j,k)[1]
                                         +intmat::Nxz(i,j,k)[0]*spins::csz(i,j,k)[0]-intmat::Nxz(i,j,k)[1]*spins::csz(i,j,k)[1];
                    fields::Hkx(i,j,k)[1]=intmat::Nxx(i,j,k)[0]*spins::csx(i,j,k)[1]+intmat::Nxx(i,j,k)[1]*spins::csx(i,j,k)[0]
                                         +intmat::Nxy(i,j,k)[0]*spins::csy(i,j,k)[1]+intmat::Nxy(i,j,k)[1]*spins::csy(i,j,k)[0]
                                         +intmat::Nxz(i,j,k)[0]*spins::csz(i,j,k)[1]+intmat::Nxz(i,j,k)[1]*spins::csz(i,j,k)[0];

                    fields::Hky(i,j,k)[0]=intmat::Nyx(i,j,k)[0]*spins::csx(i,j,k)[0]-intmat::Nyx(i,j,k)[1]*spins::csx(i,j,k)[1]
                                         +intmat::Nyy(i,j,k)[0]*spins::csy(i,j,k)[0]-intmat::Nyy(i,j,k)[1]*spins::csy(i,j,k)[1]
                                         +intmat::Nyz(i,j,k)[0]*spins::csz(i,j,k)[0]-intmat::Nyz(i,j,k)[1]*spins::csz(i,j,k)[1];
                    fields::Hky(i,j,k)[1]=intmat::Nyx(i,j,k)[0]*spins::csx(i,j,k)[1]+intmat::Nyx(i,j,k)[1]*spins::csx(i,j,k)[0]
                                         +intmat::Nyy(i,j,k)[0]*spins::csy(i,j,k)[1]+intmat::Nyy(i,j,k)[1]*spins::csy(i,j,k)[0]
                                         +intmat::Nyz(i,j,k)[0]*spins::csz(i,j,k)[1]+intmat::Nyz(i,j,k)[1]*spins::csz(i,j,k)[0];

                    fields::Hkz(i,j,k)[0]=intmat::Nzx(i,j,k)[0]*spins::csx(i,j,k)[0]-intmat::Nzx(i,j,k)[1]*spins::csx(i,j,k)[1]
                                         +intmat::Nzy(i,j,k)[0]*spins::csy(i,j,k)[0]-intmat::Nzy(i,j,k)[1]*spins::csy(i,j,k)[1]
                                         +intmat::Nzz(i,j,k)[0]*spins::csz(i,j,k)[0]-intmat::Nzz(i,j,k)[1]*spins::csz(i,j,k)[1];
                    fields::Hkz(i,j,k)[1]=intmat::Nzx(i,j,k)[0]*spins::csx(i,j,k)[1]+intmat::Nzx(i,j,k)[1]*spins::csx(i,j,k)[0]
                                         +intmat::Nzy(i,j,k)[0]*spins::csy(i,j,k)[1]+intmat::Nzy(i,j,k)[1]*spins::csy(i,j,k)[0]
                                         +intmat::Nzz(i,j,k)[0]*spins::csz(i,j,k)[1]+intmat::Nzz(i,j,k)[1]*spins::csz(i,j,k)[0];

				}
			}
		}
		fields::transformFieldsBack();
	}
}
