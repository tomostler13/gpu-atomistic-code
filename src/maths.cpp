// File: maths.cpp
// Author:Tom Ostler
// Last-modified: 22 Jan 2013 16:32:24
#include "../inc/maths.h"
#include "../inc/random.h"
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
}
