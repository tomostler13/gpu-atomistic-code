// File: fields.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 19 Aug 2022 13:12:30
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
#include "../inc/llg.h"
#ifdef CUDA
#include "../inc/cuda.h"
#endif

#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace fields
{
    Array5D<fftw_complex> Hk;
    Array5D<fftw_complex> Hr;
    Array4D<fftw_complex> hHr,hHk;
    Array4D<fftw_complex> dipHr,dipHk;
    Array<double> Hx,Hy,Hz,Hthx,Hthy,Hthz,HDemagx,HDemagy,HDemagz;
    Array<double> H4sx,H4sy,H4sz;
    Array<double> Hstagx,Hstagy,Hstagz;
    Array2D<double> sofield;
    fftw_plan HP,dHP;
    void initFields()
    {
        if(config::exchm==0)
        {
            Hk.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::Nk[2]*geom::zpdim[2]);
            Hr.resize(geom::ucm.GetNMS(),3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::Nk[2]*geom::zpdim[2]);
            //plan the transforms as in-place as we do not need to use the fft arrays
            //as we copy the data back to the normal field arrayl
            int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zps;
            int idist=geom::zps;

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of fields matrix (back transform)" << std::endl);
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
            FIXOUT(config::Log,"Direction (sign) = " << "FFTW_BACKWARD" << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            HP = fftw_plan_many_dft(3,n,geom::ucm.GetNMS()*3,Hk.ptr(),inembed,istride,idist,Hr.ptr(),onembed,ostride,odist,FFTW_BACKWARD,FFTW_ESTIMATE);
        }
        else if(config::dipm==0 && config::inc_dip==true)
        {
            dipHk.resize(3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::Nk[2]*geom::zpdim[2]);
            dipHr.resize(3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::Nk[2]*geom::zpdim[2]);
            //plan the transforms as in-place as we do not need to use the fft arrays
            //as we copy the data back to the normal field arrayl
            int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zps;
            int idist=geom::zps;

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of fields matrix (back transform)" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << 3 << std::endl);
            FIXOUT(config::Log,"Pointer of reciprocal space fields (dipHk):" << dipHk.ptr() << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUT(config::Log,"Pointer of real space fields (dipHr):" << dipHr.ptr() << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) = " << "FFTW_BACKWARD" << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            dHP = fftw_plan_many_dft(3,n,3,dipHk.ptr(),inembed,istride,idist,dipHr.ptr(),onembed,ostride,odist,FFTW_BACKWARD,FFTW_PATIENT);
        }
        if(config::exchm>98)
        {
            hHk.resize(3,geom::dim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::Nk[2]*geom::zpdim[2]);
            hHr.resize(3,geom::dim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::Nk[2]*geom::zpdim[2]);
            //plan the transforms as in-place as we do not need to use the fft arrays
            //as we copy the data back to the normal field arrayl
            int n[2]={geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zpdim[1]*geom::Nk[1]*geom::zpdim[2]*geom::Nk[2];
            int idist=geom::zpdim[1]*geom::Nk[1]*geom::zpdim[2]*geom::Nk[2];

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of fields matrix (back transform)" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[1],n[2],0);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 2 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << 3*geom::dim[0]*geom::Nk[0] << std::endl);
            FIXOUT(config::Log,"Pointer of reciprocal space fields (Hk):" << hHk.ptr() << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],0);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUT(config::Log,"Pointer of real space fields (hHr):" << hHr.ptr() << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],0);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) = " << "FFTW_BACKWARD" << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            HP = fftw_plan_many_dft(2,n,3*geom::dim[0]*geom::Nk[0],hHk.ptr(),inembed,istride,idist,hHr.ptr(),onembed,ostride,odist,FFTW_BACKWARD,FFTW_ESTIMATE);

        }
//        std::cout << "Resizing field arrays to
        Hx.resize(geom::nspins);
        Hy.resize(geom::nspins);
        Hz.resize(geom::nspins);
        HDemagx.resize(geom::nspins);
        HDemagy.resize(geom::nspins);
        HDemagz.resize(geom::nspins);
        HDemagx.IFill(0);
        HDemagy.IFill(0);
        HDemagz.IFill(0);
        H4sx.resize(geom::nspins);
        H4sy.resize(geom::nspins);
        H4sz.resize(geom::nspins);
        H4sx.IFill(0);
        H4sy.IFill(0);
        H4sz.IFill(0);
        Hstagx.resize(geom::nspins);
        Hstagy.resize(geom::nspins);
        Hstagz.resize(geom::nspins);
        Hstagx.IFill(0);
        Hstagy.IFill(0);
        Hstagz.IFill(0);
    }

    void setStagZero()
    {
        Hstagx.IFill(0);
        Hstagy.IFill(0);
        Hstagz.IFill(0);
        #ifdef CUDA
        cullg::CsetStagFieldsZero();
        #endif
    }
    void setStagField()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            unsigned int splu=geom::lu(i,3);
            fields::Hstagx(i)=fields::sofield(splu,0);
            fields::Hstagy(i)=fields::sofield(splu,1);
            fields::Hstagz(i)=fields::sofield(splu,2);
        }
        #ifdef CUDA
        cullg::CsetStagFields();
        #endif
    }
    void bfdip()
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            double h[3]={0,0,0};
            //lookup the mesh-points coords
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
                    //lookup the mesh-points coords
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
                        h[xyz]+=(3.0*sjdote*eij[xyz]*geom::mu[j]-sj[xyz]*geom::mu[j])*oomrij3;
                    }//xyz

                }
            }
            std::cout << i << "\t" << h[0]*1e-7*llg::muB << "\t" << h[1]*1e-7*llg::muB << "\t" << h[2]*1e-7*llg::muB << std::endl;
//            fields::Hx[i]=h[0]*1e-7*mat::mu*llg::muB;
//            fields::Hy[i]=h[1]*1e-7*mat::mu*llg::muB;
//            fields::Hz[i]=h[2]*1e-7*mat::mu*llg::muB;
//            std::cerr << ri[0]/geom::abc[0] << "\t" << ri[1]/geom::abc[1] << "\t" << ri[2]/geom::abc[2] << "\t" << fields::Hx[i] << "\t" << fields::Hy[i] << "\t" << fields::Hz[i] << std::endl;

        }
        exit(0);
    }//bfdip function
    //back transform the fields
    void FFTBack()
    {
        fftw_execute(HP);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            unsigned int sl=geom::sublattice[i];
            Hx[i]=Hr(sl,0,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zps));
            Hy[i]=Hr(sl,1,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zps));
            Hz[i]=Hr(sl,2,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zps));
            //std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx[i] << "\t" << fields::Hy[i] << "\t" << fields::Hz[i] << std::endl;
            //std::cin.get();
        }
    }
    void hFFTBack()
    {
        fftw_execute(HP);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Hx[i]=hHr(0,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zpdim[1]*geom::Nk[1]*geom::zpdim[2]*geom::Nk[2]));
            Hy[i]=hHr(1,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zpdim[1]*geom::Nk[1]*geom::zpdim[2]*geom::Nk[2]));
            Hz[i]=hHr(2,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zpdim[1]*geom::Nk[1]*geom::zpdim[2]*geom::Nk[2]));
            //std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx[i] << "\t" << fields::Hy[i] << "\t" << fields::Hz[i] << std::endl;
        }
        //std::cout << "*** end of fields after back transform and normalization ***" << std::endl;
        //std::cin.get();
    }
    //fourier transform method for calculating dipolar field
    void ftdip()
    {
        if(config::exchm<98)
        {
            spins::FFTForward();
            util::cpuConvFourier();
            fields::FFTBack();

        }
        else
        {
            spins::hFFTForward();
            util::hcpuConvFourier();
            fields::hFFTBack();

        }
    }
    void eftdip()
    {
        if(config::exchm<98)
        {
            spins::eFFTForward();
            util::cpuConvFourier();
            fields::FFTBack();
        }
        else
        {
            spins::heFFTForward();
            util::hcpuConvFourier();
            fields::hFFTBack();
        }
    }
    void dipFFTBack()
    {
        fftw_execute(dHP);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            HDemagx[i]=dipHr(0,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zps));
            HDemagy[i]=dipHr(1,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zps));
            HDemagz[i]=dipHr(2,lc[0],lc[1],lc[2])[0]/(static_cast<double>(geom::zps));
            //std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx[i] << "\t" << fields::Hy[i] << "\t" << fields::Hz[i] << std::endl;
            //std::cin.get();
        }
    }
    //fourier transform method for calculating dipolar field
    void dipftdip()
    {
        spins::dipFFTForward();
        util::dipcpuConvFourier();
        fields::dipFFTBack();
    }
    void dipeftdip()
    {
        spins::dipeFFTForward();
        util::dipcpuConvFourier();
        fields::dipFFTBack();
    }

}
