// File: fields.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 26 Sep 2014 13:49:22
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
#include "../inc/spins.h"
#include "../inc/util.h"
#include "../inc/intmat.h"
#include "../inc/defines.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace fields
{
    Array5D<fftw_complex> Hk;
    Array5D<double> Hr;
    Array<double> Hx,Hy,Hz,Hthx,Hthy,Hthz;
    fftw_plan HP;
    void initFields(int argc,char *argv[])
    {
        Hk.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Hr.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::Nk[2]*geom::zpdim[2]);
        Hx.resize(geom::nspins);
        Hy.resize(geom::nspins);
        Hz.resize(geom::nspins);
        //plan the transforms as in-place as we do not need to use the fft arrays
        //as we copy the data back to the normal field arrayl
        int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
        int *inembed=n;
        int *onembed=n;
        int istride=1;
        int ostride=1;
        int odist=geom::zps;
        int idist=geom::cplxdim;

        config::openLogFile();
        config::printline(config::Log);
        FIXOUT(config::Log,"Parameters entering into c2r FFTW plan of fields matrix (back transform)" << std::endl);
        FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
        FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
        FIXOUT(config::Log,"How many (FFT's) = " << geom::ucm.GetNMS()*3 << std::endl);
        FIXOUT(config::Log,"Pointer of reciprocal space fields (Hk):" << Hk.ptr() << std::endl);
        FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
        FIXOUT(config::Log,"istride = " << istride << std::endl);
        FIXOUT(config::Log,"idist = " << idist << std::endl);
        FIXOUT(config::Log,"Pointer of real space fields (Hr):" << Hr.ptr() << std::endl);
        FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
        FIXOUT(config::Log,"ostride = " << ostride << std::endl);
        FIXOUT(config::Log,"odist = " << odist << std::endl);
        FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
        HP = fftw_plan_many_dft_c2r(3,n,geom::ucm.GetNMS()*3,Hk.ptr(),inembed,istride,idist,Hr.ptr(),onembed,ostride,odist,FFTW_PATIENT);
    }

    /*void bfdip()
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
//            fields::Hx[i]=h[0]*1e-7*mat::mu*mat::muB;
//            fields::Hy[i]=h[1]*1e-7*mat::mu*mat::muB;
//            fields::Hz[i]=h[2]*1e-7*mat::mu*mat::muB;
            std::cerr << ri[0]/geom::abc[0] << "\t" << ri[1]/geom::abc[1] << "\t" << ri[2]/geom::abc[2] << "\t" << fields::Hx[i] << "\t" << fields::Hy[i] << "\t" << fields::Hz[i] << std::endl;

        }
    }//bfdip function*/
    //back transform the fields
    void FFTBack()
    {
        fftw_execute(HP);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            unsigned int sl=geom::sublattice[i];
            Hx[i]=Hr(sl,0,lc[0],lc[1],lc[2])/(static_cast<double>(geom::zps));
            Hy[i]=Hr(sl,1,lc[0],lc[1],lc[2])/(static_cast<double>(geom::zps));
            Hz[i]=Hr(sl,2,lc[0],lc[1],lc[2])/(static_cast<double>(geom::zps));
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
