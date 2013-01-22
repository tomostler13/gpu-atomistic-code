// File: fields.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 22 Jan 2013 10:57:54
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
#include "../inc/spins.h"
#include "../inc/util.h"
#include "../inc/intmat.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace fields
{
    Array3D<fftw_complex> Hkx,Hky,Hkz;
    Array3D<double> Hrx,Hry,Hrz;
    Array<double> Hx,Hy,Hz,Hthx,Hthy,Hthz;
    fftw_plan HxP,HyP,HzP;
    void initFields(int argc,char *argv[])
    {
        Hkx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Hky.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Hkz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Hx.resize(geom::nspins);
        Hy.resize(geom::nspins);
        Hz.resize(geom::nspins);
        Hrx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Hry.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Hrz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        //plan the transforms as in-place as we do not need to use the fft arrays
        //as we copy the data back to the normal field arrayl
        HxP = fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Hkx.ptr(),Hrx.ptr(),FFTW_ESTIMATE);
        HyP = fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Hky.ptr(),Hry.ptr(),FFTW_ESTIMATE);
        HzP = fftw_plan_dft_c2r_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Hkz.ptr(),Hrz.ptr(),FFTW_ESTIMATE);
    }

    void bfdip()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            double h[3]={0,0,0};
            //lookup the k-points coords
            unsigned int lci[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            double ri[3]={0,0,0};
            //spin components
            double si[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            //determine the real space vector from origin to atom i
            for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
            {
                ri[xyz]=geom::abc[xyz]*double(lci[xyz])/double(geom::Nk[xyz]);
            }
            for(unsigned int j = 0 ; j < geom::nspins ; j++)
            {
                if(i!=j)
                {
                    //lookup the k-points coords
                    unsigned int lcj[3]={geom::lu(j,0),geom::lu(j,1),geom::lu(j,2)};
                    double rj[3]={0,0,0};
                    double sj[3]={spins::Sx[j],spins::Sy[j],spins::Sz[j]};
                    double rij[3]={0,0,0};
                    //determine the real space vector from origin to atom i
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        rj[xyz]=geom::abc[xyz]*double(lcj[xyz])/double(geom::Nk[xyz]);
                        rij[xyz]=rj[xyz]-ri[xyz];
                    }

                    double mrij=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                    double oomrij3=1./(mrij*mrij*mrij);
                    double eij[3]={rij[0]/mrij,rij[1]/mrij,rij[2]/mrij};
                    double sjdote=sj[0]*eij[0]+sj[1]*eij[1]+sj[2]*eij[2];


                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        h[xyz]+=(3.0*sjdote*eij[xyz]-sj[xyz])*oomrij3;
                    }//xyz

                }
            }
            fields::Hx[i]=h[0]*1e-7*mat::mu*mat::muB;
            fields::Hy[i]=h[1]*1e-7*mat::mu*mat::muB;
            fields::Hz[i]=h[2]*1e-7*mat::mu*mat::muB;
            std::cerr << ri[0]/geom::abc[0] << "\t" << ri[1]/geom::abc[1] << "\t" << ri[2]/geom::abc[2] << "\t" << fields::Hx[i] << "\t" << fields::Hy[i] << "\t" << fields::Hz[i] << std::endl;

        }
    }//bfdip function
    //back transform the fields
    void FFTBack()
    {
        fftw_execute(HxP);
        fftw_execute(HyP);
        fftw_execute(HzP);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Hx[i]=Hrx(lc[0],lc[1],lc[2])/(double(geom::zps));
            Hy[i]=Hry(lc[0],lc[1],lc[2])/(double(geom::zps));
            Hz[i]=Hrz(lc[0],lc[1],lc[2])/(double(geom::zps));
            //std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx[i] << "\t" << fields::Hy[i] << "\t" << fields::Hz[i] << std::endl;
            //std::cin.get();
        }
    }
    //fourier transform method for calculating dipolar field
    void ftdip()
    {
        spins::FFTForward();
        util::cpuConvFourier();
        fields::FFTBack();
    }
    void eftdip()
    {
        spins::eFFTForward();
        util::cpuConvFourier();
        fields::FFTBack();
    }
}
