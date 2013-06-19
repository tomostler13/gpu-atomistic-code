// File: llg.cpp
// Author:Tom Ostler
// Created: 21 Jan 2013
// Last-modified: 19 Jun 2013 17:51:31
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include <omp.h>
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/intmat.h"
#include "../inc/defines.h"
#include "../inc/llg.h"
#include "../inc/random.h"
#include "../inc/mat.h"
#include "../inc/exch.h"
#include "../inc/anis.h"
namespace llgCPU
{
    Array<double> fnx;
    Array<double> fny;
    Array<double> fnz;
    void initLLG(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**LLG CPU details***" << std::endl;
        FIXOUT(config::Info,"Initializing arrays:" << std::flush);
        fnx.resize(geom::nspins);
        fny.resize(geom::nspins);
        fnz.resize(geom::nspins);
        fnx.IFill(0);
        fny.IFill(0);
        fnz.IFill(0);

        SUCCESS(config::Info);
        FIXOUT(config::Info,"Resizing llg (CPU) work arrays:" << std::flush);
        spins::eSx.resize(geom::nspins);
        spins::eSy.resize(geom::nspins);
        spins::eSz.resize(geom::nspins);
        fields::Hthx.resize(geom::nspins);
        fields::Hthy.resize(geom::nspins);
        fields::Hthz.resize(geom::nspins);
        SUCCESS(config::Info);
    }
    //for uniform temperature with interaction matrix for field calculation
    void llgCPU(unsigned int t,bool& faf)
    {
        //calculate the 2 spin fields (dipolar, exchange, anisotropy)
        fields::ftdip();
        if(llg::rscf==true)
        {
            //spins::calcRealSpaceCorrelationFunction();
        }
/*        //FOR DEBUGGING THE FIELD
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
        }
        exit(0);*/

        const double sqrtT=sqrt(llg::T);
        #pragma omp parallel for private (i) shared(fields::Hthx,fields::Hthy,fields::Hthz,sqrtT,mat::sigma)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            fields::Hthx[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthy[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthz[i]=sqrtT*mat::sigma[i]*Random::normal();

        }
        #pragma omp parallel for private (i) shared(spins::eSx,spins::eSy,spins::eSz,fnx,fny,fnz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            fnx[i]=llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]);
            fny[i]=llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]);
            fnz[i]=llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2]);
            spins::eSx[i] = s[0] + fnx[i]*llg::rdt;
            spins::eSy[i] = s[1] + fny[i]*llg::rdt;
            spins::eSz[i] = s[2] + fnz[i]*llg::rdt;
			const double mods = sqrt(spins::eSx[i]*spins::eSx[i]+spins::eSy[i]*spins::eSy[i]+spins::eSz[i]*spins::eSz[i]);
			spins::eSx[i]/=mods;
			spins::eSy[i]/=mods;
			spins::eSz[i]/=mods;
        }
        //perform the calculation of the 2 spin fields using the euler spin arrays
        fields::eftdip();
        #pragma omp parallel for private (i) shared(spins::Sx,spins::Sy,spins::Sz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            const double fnp1[3]={llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]),llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]),llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2])};

            spins::Sx[i] +=(0.5*(fnx[i]+fnp1[0])*llg::rdt);
            spins::Sy[i] +=(0.5*(fny[i]+fnp1[1])*llg::rdt);
            spins::Sz[i] +=(0.5*(fnz[i]+fnp1[2])*llg::rdt);
			const double mods = sqrt(spins::Sx[i]*spins::Sx[i]+spins::Sy[i]*spins::Sy[i]+spins::Sz[i]*spins::Sz[i]);
			spins::Sx[i]/=mods;
			spins::Sy[i]/=mods;
			spins::Sz[i]/=mods;
        }
    }
    //for uniform temperature with interaction matrix for field calculation (uniform field)
    void llgCPU(unsigned int t)
    {
        //calculate the 2 spin fields (dipolar, exchange, anisotropy)
        fields::ftdip();
        if(llg::rscf==true)
        {
            //spins::calcRealSpaceCorrelationFunction();
        }
/*        //FOR DEBUGGING THE FIELD
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
        }
        exit(0);*/

        const double sqrtT=sqrt(llg::T);
        #pragma omp parallel for private (i) shared(fields::Hthx,fields::Hthy,fields::Hthz,sqrtT,mat::sigma)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            fields::Hthx[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthy[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthz[i]=sqrtT*mat::sigma[i]*Random::normal();

        }
        #pragma omp parallel for private (i) shared(spins::eSx,spins::eSy,spins::eSz,fnx,fny,fnz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            double h[3]={llg::applied[0]+fields::Hthx[i]+fields::Hx[i],llg::applied[1]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+llg::applied[2]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            fnx[i]=llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]);
            fny[i]=llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]);
            fnz[i]=llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2]);
            spins::eSx[i] = s[0] + fnx[i]*llg::rdt;
            spins::eSy[i] = s[1] + fny[i]*llg::rdt;
            spins::eSz[i] = s[2] + fnz[i]*llg::rdt;
			const double mods = sqrt(spins::eSx[i]*spins::eSx[i]+spins::eSy[i]*spins::eSy[i]+spins::eSz[i]*spins::eSz[i]);
			spins::eSx[i]/=mods;
			spins::eSy[i]/=mods;
			spins::eSz[i]/=mods;
        }
        //perform the calculation of the 2 spin fields using the euler spin arrays
        fields::eftdip();
        #pragma omp parallel for private (i) shared(spins::Sx,spins::Sy,spins::Sz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={llg::applied[0]+fields::Hthx[i]+fields::Hx[i],llg::applied[1]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+llg::applied[2]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            const double fnp1[3]={llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]),llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]),llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2])};

            spins::Sx[i] +=(0.5*(fnx[i]+fnp1[0])*llg::rdt);
            spins::Sy[i] +=(0.5*(fny[i]+fnp1[1])*llg::rdt);
            spins::Sz[i] +=(0.5*(fnz[i]+fnp1[2])*llg::rdt);
			const double mods = sqrt(spins::Sx[i]*spins::Sx[i]+spins::Sy[i]*spins::Sy[i]+spins::Sz[i]*spins::Sz[i]);
			spins::Sx[i]/=mods;
			spins::Sy[i]/=mods;
			spins::Sz[i]/=mods;
        }
    }
    //for on-site temperature and using interaction matrix for field calculation
    void llgCPU(unsigned int t,Array<double>& T)
    {
        //calculate the 2 spin fields (dipolar, exchange, anisotropy)
        fields::ftdip();
        if(llg::rscf==true)
        {
            //spins::calcRealSpaceCorrelationFunction();
        }
        //FOR DEBUGGING THE FIELD
/*        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
        }
        exit(0);*/

        #pragma omp parallel for private (i) shared(fields::Hthx,fields::Hthy,fields::Hthz,sqrtT,mat::sigma)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double sqrtT=sqrt(T[i]);
            fields::Hthx[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthy[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthz[i]=sqrtT*mat::sigma[i]*Random::normal();

        }
        #pragma omp parallel for private (i) shared(spins::eSx,spins::eSy,spins::eSz,fnx,fny,fnz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            fnx[i]=llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]);
            fny[i]=llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]);
            fnz[i]=llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2]);
            spins::eSx[i] = s[0] + fnx[i]*llg::rdt;
            spins::eSy[i] = s[1] + fny[i]*llg::rdt;
            spins::eSz[i] = s[2] + fnz[i]*llg::rdt;
			const double mods = sqrt(spins::eSx[i]*spins::eSx[i]+spins::eSy[i]*spins::eSy[i]+spins::eSz[i]*spins::eSz[i]);
			spins::eSx[i]/=mods;
			spins::eSy[i]/=mods;
			spins::eSz[i]/=mods;
        }
        //perform the calculation of the 2 spin fields using the euler spin arrays
        fields::eftdip();
        #pragma omp parallel for private (i) shared(spins::Sx,spins::Sy,spins::Sz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            const double fnp1[3]={llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]),llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]),llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2])};

            spins::Sx[i] +=(0.5*(fnx[i]+fnp1[0])*llg::rdt);
            spins::Sy[i] +=(0.5*(fny[i]+fnp1[1])*llg::rdt);
            spins::Sz[i] +=(0.5*(fnz[i]+fnp1[2])*llg::rdt);
			const double mods = sqrt(spins::Sx[i]*spins::Sx[i]+spins::Sy[i]*spins::Sy[i]+spins::Sz[i]*spins::Sz[i]);
			spins::Sx[i]/=mods;
			spins::Sy[i]/=mods;
			spins::Sz[i]/=mods;
        }
    }
    //Exchange field calculated using a CSR neighbbour list with on-site temperature
    void llgCPU(unsigned int t,Array<double>& T,Array<unsigned int>& xadj,Array<unsigned int>& adjncy)
    {
        //calculate the 2 spin fields (dipolar, exchange, anisotropy)
        if(config::incdip)
        {
            fields::bfdip();
        }
        if(llg::rscf==true)
        {
            //spins::calcRealSpaceCorrelationFunction();
        }
        //FOR DEBUGGING THE FIELD
/*        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
        }
        exit(0);*/

        #pragma omp parallel for private (i) shared(fields::Hthx,fields::Hthy,fields::Hthz,sqrtT,mat::sigma)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double sqrtT=sqrt(T[i]);
            fields::Hthx[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthy[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthz[i]=sqrtT*mat::sigma[i]*Random::normal();
            //std::cout << fields::Hthx[i] << "\t" << fields::Hthy[i] << "\t" << fields::Hthz[i] << std::endl;
        }
        //std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        #pragma omp parallel for private (i) shared(spins::eSx,spins::eSy,spins::eSz,fnx,fny,fnz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int spec=mat::speclist[i];
            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};

            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};

            for(unsigned int n=xadj[i] ; n < xadj[i+1] ; n++)
            {
                unsigned int neigh=adjncy[n];
                h[0]+=exch::Jxx[neigh]*spins::Sx[neigh];
                h[1]+=exch::Jyy[neigh]*spins::Sy[neigh];
                h[2]+=exch::Jzz[neigh]*spins::Sz[neigh];
            }

            const double sdotn=s[0]*anis::uniaxial_unit(spec,0)+s[1]*anis::uniaxial_unit(spec,1)+s[2]*anis::uniaxial_unit(spec,2);
            h[0]+=anis::dT(spec,0,0)*sdotn;
            h[1]+=anis::dT(spec,1,1)*sdotn;
            h[2]+=anis::dT(spec,2,2)*sdotn;
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};
            fnx[i]=llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]);
            fny[i]=llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]);
            fnz[i]=llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2]);
            spins::eSx[i] = s[0] + fnx[i]*llg::rdt;
            spins::eSy[i] = s[1] + fny[i]*llg::rdt;
            spins::eSz[i] = s[2] + fnz[i]*llg::rdt;


			const double mods = sqrt(spins::eSx[i]*spins::eSx[i]+spins::eSy[i]*spins::eSy[i]+spins::eSz[i]*spins::eSz[i]);
			spins::eSx[i]/=mods;
			spins::eSy[i]/=mods;
			spins::eSz[i]/=mods;
        }
        #pragma omp parallel for private (i) shared(spins::Sx,spins::Sy,spins::Sz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int spec=mat::speclist[i];
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};
            for(unsigned int n=xadj[i] ; n < xadj[i+1] ; n++)
            {
                unsigned int neigh=adjncy[n];
                h[0]+=exch::Jxx[neigh]*spins::eSx[neigh];
                h[1]+=exch::Jyy[neigh]*spins::eSy[neigh];
                h[2]+=exch::Jzz[neigh]*spins::eSz[neigh];
            }

            const double sdotn=s[0]*anis::uniaxial_unit(spec,0)+s[1]*anis::uniaxial_unit(spec,1)+s[2]*anis::uniaxial_unit(spec,2);
            h[0]+=anis::dT(spec,0,0)*sdotn;
            h[1]+=anis::dT(spec,1,1)*sdotn;
            h[2]+=anis::dT(spec,2,2)*sdotn;

            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            const double fnp1[3]={llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]),llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]),llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2])};

            spins::Sx[i] +=(0.5*(fnx[i]+fnp1[0])*llg::rdt);
            spins::Sy[i] +=(0.5*(fny[i]+fnp1[1])*llg::rdt);
            spins::Sz[i] +=(0.5*(fnz[i]+fnp1[2])*llg::rdt);
			const double mods = sqrt(spins::Sx[i]*spins::Sx[i]+spins::Sy[i]*spins::Sy[i]+spins::Sz[i]*spins::Sz[i]);
			spins::Sx[i]/=mods;
			spins::Sy[i]/=mods;
			spins::Sz[i]/=mods;
        }
    }
    //Exchange field calculated using a CSR neighbbour list with uniform temperature
    void llgCPU(unsigned int t,Array<unsigned int>& xadj,Array<unsigned int>& adjncy)
    {
        //calculate the 2 spin fields (dipolar, exchange, anisotropy)
        if(config::incdip)
        {
            fields::bfdip();
        }
        if(llg::rscf==true)
        {
            //spins::calcRealSpaceCorrelationFunction();
        }
        //FOR DEBUGGING THE FIELD
/*        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
        }
        exit(0);*/

        #pragma omp parallel for private (i) shared(fields::Hthx,fields::Hthy,fields::Hthz,sqrtT,mat::sigma)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

            const double sqrtT=llg::T;
            fields::Hthx[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthy[i]=sqrtT*mat::sigma[i]*Random::normal();
            fields::Hthz[i]=sqrtT*mat::sigma[i]*Random::normal();

        }
        #pragma omp parallel for private (i) shared(spins::eSx,spins::eSy,spins::eSz,fnx,fny,fnz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const unsigned int spec=mat::speclist[i];
            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};
            for(unsigned int n=xadj[i] ; n < xadj[i+1] ; n++)
            {
                unsigned int neigh=adjncy[n];
                h[0]+=exch::Jxx[neigh]*spins::Sx[neigh];
                h[1]+=exch::Jyy[neigh]*spins::Sy[neigh];
                h[2]+=exch::Jzz[neigh]*spins::Sz[neigh];
            }
            const double sdotn=s[0]*anis::uniaxial_unit(spec,0)+s[1]*anis::uniaxial_unit(spec,1)+s[2]*anis::uniaxial_unit(spec,2);
            h[0]+=anis::dT(spec,0,0)*sdotn;
            h[1]+=anis::dT(spec,1,1)*sdotn;
            h[2]+=anis::dT(spec,2,2)*sdotn;
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            fnx[i]=llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]);
            fny[i]=llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]);
            fnz[i]=llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2]);
            spins::eSx[i] = s[0] + fnx[i]*llg::rdt;
            spins::eSy[i] = s[1] + fny[i]*llg::rdt;
            spins::eSz[i] = s[2] + fnz[i]*llg::rdt;
			const double mods = sqrt(spins::eSx[i]*spins::eSx[i]+spins::eSy[i]*spins::eSy[i]+spins::eSz[i]*spins::eSz[i]);
			spins::eSx[i]/=mods;
			spins::eSy[i]/=mods;
			spins::eSz[i]/=mods;
        }
        #pragma omp parallel for private (i) shared(spins::Sx,spins::Sy,spins::Sz)
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const unsigned int spec=mat::speclist[i];
            const double s[3]={spins::eSx[i],spins::eSy[i],spins::eSz[i]};
            double h[3]={fields::HAppx[i]+fields::Hthx[i]+fields::Hx[i],fields::HAppy[i]+fields::Hthy[i]+fields::Hy[i],fields::Hz[i]+fields::Hthz[i]+fields::HAppz[i]};
            for(unsigned int n=xadj[i] ; n < xadj[i+1] ; n++)
            {
                unsigned int neigh=adjncy[n];
                h[0]+=exch::Jxx[neigh]*spins::eSx[neigh];
                h[1]+=exch::Jyy[neigh]*spins::eSy[neigh];
                h[2]+=exch::Jzz[neigh]*spins::eSz[neigh];
            }
            const double sdotn=s[0]*anis::uniaxial_unit(spec,0)+s[1]*anis::uniaxial_unit(spec,1)+s[2]*anis::uniaxial_unit(spec,2);
            h[0]+=anis::dT(spec,0,0)*sdotn;
            h[1]+=anis::dT(spec,1,1)*sdotn;
            h[2]+=anis::dT(spec,2,2)*sdotn;
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};

            const double fnp1[3]={llg::llgpf[i]*(sxh[0]+mat::lambda[i]*sxsxh[0]),llg::llgpf[i]*(sxh[1]+mat::lambda[i]*sxsxh[1]),llg::llgpf[i]*(sxh[2]+mat::lambda[i]*sxsxh[2])};

            spins::Sx[i] +=(0.5*(fnx[i]+fnp1[0])*llg::rdt);
            spins::Sy[i] +=(0.5*(fny[i]+fnp1[1])*llg::rdt);
            spins::Sz[i] +=(0.5*(fnz[i]+fnp1[2])*llg::rdt);
			const double mods = sqrt(spins::Sx[i]*spins::Sx[i]+spins::Sy[i]*spins::Sy[i]+spins::Sz[i]*spins::Sz[i]);
			spins::Sx[i]/=mods;
			spins::Sy[i]/=mods;
			spins::Sz[i]/=mods;
        }
    }


}
