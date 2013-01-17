// File: fields.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 17 Jan 2013 14:53:25
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
    Array<double> Hx,Hy,Hz;
    fftw_plan HxP,HyP,HzP;
    void initFields(int argc,char *argv[])
    {
        Hkx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Hky.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Hkz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Hx.resize(geom::nspins);
        Hy.resize(geom::nspins);
        Hz.resize(geom::nspins);
        //plan the transforms as in-place as we do not need to use the fft arrays
        //as we copy the data back to the normal field arrayl
        HxP = fftw_plan_dft_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Hkx.ptr(),Hkx.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
        HyP = fftw_plan_dft_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Hky.ptr(),Hky.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
        HzP = fftw_plan_dft_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Hkz.ptr(),Hkz.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
    }

    void bfdip()
    {
        for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                {
                    for(unsigned int t = 0 ; t < geom::nauc ; t++)
                    {
                        double h[3]={0,0,0};
                        double ri[3]={0,0,0};
                        double qi[3]={double(i),double(j),double(k)};
                        for(unsigned int a = 0 ; a < 3 ; a++)
                        {
                            for(unsigned int b = 0 ; b < 3 ; b++)
                            {
                                ri[a]+=geom::L(a,b)*qi[b];
                            }
                            ri[a]+=geom::ucm.GetComp(t,a);
                            ri[a]*=geom::abc[a];
                        }
                        //lookup the spin number
                        unsigned int sni=(i*geom::dim[2]*geom::dim[1]+j*geom::dim[2]+k)*geom::nauc+t;
                        double si[3]={spins::Sx[sni],spins::Sy[sni],spins::Sz[sni]};

                        for(unsigned int x = 0 ; x < geom::dim[0] ; x++)
                        {
                            for(unsigned int y = 0 ; y < geom::dim[1] ; y++)
                            {
                                for(unsigned int z = 0 ; z < geom::dim[2] ; z++)
                                {
                                    for(unsigned int t2 = 0 ; t2 < geom::nauc ; t2++)
                                    {
                                        if(i==x && j==y && k==z && t==t2)
                                        {
                                            h[0]+=0;
                                            h[1]+=0;
                                            h[2]+=0;
                                        }
                                        else
                                        {
                                            double rj[3]={0,0,0};
                                            double qj[3]={double(x),double(y),double(z)};
                                            for(unsigned int c = 0 ; c < 3 ; c++)
                                            {
                                                for(unsigned int d = 0 ; d < 3 ; d++)
                                                {
                                                    rj[c]+=geom::L(c,d)*qj[d];
                                                }//d
                                                rj[c]+=geom::ucm.GetComp(t2,c);
                                                rj[c]*=geom::abc[c];
                                            }//d
                                            //lookup the spin number
                                            unsigned int snj=(x*geom::dim[2]*geom::dim[1]+y*geom::dim[2]+z)*geom::nauc+t2;
                                            double sj[3]={spins::Sx[snj],spins::Sy[snj],spins::Sz[snj]};
                                            double rij[3]={0,0,0};
                                            for(unsigned int comp = 0 ; comp < 3 ; comp++)
                                            {
                                                rij[comp]=rj[comp]-ri[comp];
                                            }//comp
                                            double mrij=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                                            double oomrij3=1./(mrij*mrij*mrij);
                                            double eij[3]={rij[0]/mrij,rij[1]/mrij,rij[2]/mrij};
                                            double sjdote=sj[0]*eij[0]+sj[1]*eij[1]+sj[2]*eij[2];


                                            for(unsigned int comp = 0 ; comp < 3 ; comp++)
                                            {
                                                h[comp]+=(3.0*sjdote*eij[comp]-sj[comp])*oomrij3;
                                            }//comp

                                        }
                                    }//t2
                                }//z
                            }//y
                        }//x
                        fields::Hx[sni]=h[0]*1e-7*mat::mu*mat::muB;
                        fields::Hy[sni]=h[1]*1e-7*mat::mu*mat::muB;
                        fields::Hz[sni]=h[2]*1e-7*mat::mu*mat::muB;
                        std::cout << ri[0]/geom::abc[0] << "\t" << ri[1]/geom::abc[1] << "\t" << ri[2]/geom::abc[2] << "\t" << fields::Hx[sni] << "\t" << fields::Hy[sni] << "\t" << fields::Hz[sni] << std::endl;

                    }//t
                }//k
            }//j
        }//i
    }//bfdip function
    //back transform the fields
    void FFTBack()
    {
        fftw_execute(HxP);
        fftw_execute(HyP);
        fftw_execute(HzP);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            Hx[i]=Hkx.ReturnDirectReal(intmat::zpsn[i])/double(geom::zps);
            Hy[i]=Hky.ReturnDirectReal(intmat::zpsn[i])/double(geom::zps);
            Hz[i]=Hkz.ReturnDirectReal(intmat::zpsn[i])/double(geom::zps);
        }
    }
    //fourier transform method for calculating dipolar field
    void ftdip()
    {
        spins::FFTForward();
        util::cpuConvFourier();
        fields::FFTBack();
        unsigned int atom_counter=0;
        for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                {
                    for(unsigned int t = 0 ; t < geom::nauc ; t++)
                    {
                        double ri[3]={0,0,0};
                        double qi[3]={double(i),double(j),double(k)};
                        for(unsigned int a = 0 ; a < 3 ; a++)
                        {
                            for(unsigned int b = 0 ; b < 3 ; b++)
                            {
                                ri[a]+=geom::L(a,b)*qi[b];
                            }
                            ri[a]+=geom::ucm.GetComp(t,a);
                            ri[a]*=geom::abc[a];
                        }
                        std::cout << ri[0]/geom::abc[0] << "\t" << ri[1]/geom::abc[1] << "\t" << ri[2]/geom::abc[2] << "\t" << Hx[atom_counter] << "\t" << Hy[atom_counter] << "\t" << Hz[atom_counter] << std::endl;
                        atom_counter++;
                    }
                }
            }
        }

    }
}
