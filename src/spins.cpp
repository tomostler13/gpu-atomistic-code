// File: spins.cpp
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 25 Sep 2014 14:40:06
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
#include "../inc/intmat.h"
#include "../inc/defines.h"
#include "../inc/maths.h"
namespace spins
{
    Array5D<fftw_complex> Sk;
    Array5D<double> Sr;
    Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
    fftw_plan SP;
    unsigned int update=0;
    std::ifstream sfs;
    void initSpins(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Spin details***" << std::endl;
        FIXOUT(config::Info,"Resizing arrays:" << std::flush);
        Sk.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Sr.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Sr.IFill(0);
        Sk.IFill(0);
        Sx.resize(geom::nspins);
        Sy.resize(geom::nspins);
        Sz.resize(geom::nspins);
        Sx.IFill(geom::nspins);
        Sy.IFill(geom::nspins);
        Sz.IFill(geom::nspins);
        SUCCESS(config::Info);

        FIXOUT(config::Info,"Planning r2c and c2r transforms:" << std::flush);
        int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
        int *inembed=n;
        int *onembed=n;
        int istride=1;
        int ostride=1;
        int odist=geom::cplxdim;
        int idist=geom::zps;

        config::openLogFile();
        config::printline(config::Log);
        FIXOUT(config::Log,"Parameters entering into r2c FFTW plan of spins (forward transform)" << std::endl);
        FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
        FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
        FIXOUT(config::Log,"How many (FFT's) = " << geom::ucm.GetNMS()*3 << std::endl);
        FIXOUT(config::Log,"Pointer of real space fields (Sr):" << Sr.ptr() << std::endl);
        FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
        FIXOUT(config::Log,"istride = " << istride << std::endl);
        FIXOUT(config::Log,"idist = " << idist << std::endl);
        FIXOUT(config::Log,"Pointer of reciprocal space fields (Sk):" << Sk.ptr() << std::endl);
        FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
        FIXOUT(config::Log,"ostride = " << ostride << std::endl);
        FIXOUT(config::Log,"odist = " << odist << std::endl);
        FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
        SP = fftw_plan_many_dft_r2c(3,n,geom::ucm.GetNMS()*3,Sr.ptr(),inembed,istride,idist,Sk.ptr(),onembed,ostride,odist,FFTW_PATIENT);
        //forward transform of spin arrays
        SUCCESS(config::Info);

    }
    void FFTForward()
    {
        Sk.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            unsigned int sl=geom::sublattice[i];
            Sr(sl,0,lc[0],lc[1],lc[2])=Sx[i];
            Sr(sl,1,lc[0],lc[1],lc[2])=Sx[i];
            Sr(sl,2,lc[0],lc[1],lc[2])=Sx[i];
            //            std::cout << "Lookup\t" << lc[0] << "\t" << lc[1] << "\t" << lc[2] << "\t" << Skx(lc[0],lc[1],lc[2])[0] << "\t" << Sky(lc[0],lc[1],lc[2])[0] << "\t" << Skz(lc[0],lc[1],lc[2])[0] << std::endl;
        }
        fftw_execute(SP);
    }
    void eFFTForward()
    {
        Sk.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            unsigned int sl=geom::sublattice[i];
            Sr(sl,0,lc[0],lc[1],lc[2])=eSx[i];
            Sr(sl,1,lc[0],lc[1],lc[2])=eSy[i];
            Sr(sl,2,lc[0],lc[1],lc[2])=eSz[i];
            //            std::cout << "Lookup\t" << lc[0] << "\t" << lc[1] << "\t" << lc[2] << "\t" << Skx(lc[0],lc[1],lc[2])[0] << "\t" << Sky(lc[0],lc[1],lc[2])[0] << "\t" << Skz(lc[0],lc[1],lc[2])[0] << std::endl;
        }
        fftw_execute(SP);
    }
}
