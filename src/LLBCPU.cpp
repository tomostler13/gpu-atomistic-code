// File: LLBCPU.cpp
// Author:Tom Ostler
// Last-modified: 10 Jan 2013 13:49:54
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include "../inc/LLBCPU.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/maths.h"
#include "../inc/neigh.h"
#include "../inc/array2d.h"
#include "../inc/intmat.h"
#include "../inc/util.h"
#include "../inc/tdp.h"
#include "../inc/mat.h"
#include "../inc/sim.h"
#include "../inc/random.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace llb
{
    Array<double> fnx;
    Array<double> fny;
    Array<double> fnz;
    bool LLBi=false;

    void initLLBCPU(int argc,char *argv[])
    {
        assert(geom::gi);
        config::printline(config::Info);

		config::Info.width(45);config::Info << std::right << "*" << "**LLB CPU details***" << std::endl;

        FIXOUT(config::Info,"Initializing arrays:" << std::flush);
        fnx.resize(geom::ss);
        fny.resize(geom::ss);
        fnz.resize(geom::ss);
        fnx.IFill(0);
        fny.IFill(0);
        fnz.IFill(0);

        config::Info << "Done" << std::endl;
        LLBi=true;
    }

    void LLBCPU(unsigned int& t)
    {
        /*      FOR DEBUGGING THE NEIGHBOUR LIST
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            std::cout << "Atom i = " << i << " interacts with" << std::endl;
            for(unsigned int j = neigh::xadj[i] ; j < neigh::xadj[i+1] ; j++)
            {
                std::cout << neigh::adjncy[j] << std::endl;
            }
        }*/

        assert(LLBi);
        assert(geom::gi);
        assert(config::lcf);
        assert(neigh::ni);
        assert(util::ui);
        assert(spins::si);
        assert(fields::fi);
        assert(tdp::tdpi);
        //perform the calculation of the dipolar field
        //all methods should return the field to the
        //real part of the Hr{x,y,z} arrays. These arrays
        //should be initialized to zero before the
        //calculation begins.
        if(fields::checkdipolar==true)
        {
            if(t==0 || t%fields::dfu==0)
            {
                if(fields::dfc=="fft")
                {
                    maths::cpuConvFourier();
                }
                else if(fields::dfc=="bruteforce")
                {
                    fields::bfdip();
                }
            }
        }
        else
        {
            for(unsigned int i = 0 ; i < geom::ss ; i++)
            {
                fields::Hx[i]=0;
                fields::Hy[i]=0;
                fields::Hz[i]=0;
            }
        }
        //FOR DEBUGGING THE DIPOLAR FIELD
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            std::cout << geom::lu(i,0) << "\t" << geom::lu(i,1) << "\t" << geom::lu(i,2) << "\t" << fields::Hx(i) << "\t" << fields::Hy(i) << "\t" << fields::Hz(i) << std::endl;
        }
        exit(0);
        for(unsigned int t2=0;t2<fields::dfu;t2++)
        {
            //set the system temperature
            tdp::Tcalc();
            tdp::chiperpfp();
            tdp::chiparfp();
            tdp::exchstifffp();
            tdp::mefp();
            tdp::calcalphas();

            if(fields::checkthermal==true)
            {
                for(unsigned int i = 0 ; i < geom::ss ; i++)
                {

                    fields::GW1x[i]=Random::normal();
                    fields::GW1y[i]=Random::normal();
                    fields::GW1z[i]=Random::normal();
                    fields::GW2x[i]=Random::normal();
                    fields::GW2y[i]=Random::normal();
                    fields::GW2z[i]=Random::normal();
                }
            }
            else
            {
                for(unsigned int i = 0 ; i < geom::ss ; i++)
                {
                    fields::GW1x[i]=0.0;
                    fields::GW1y[i]=0.0;
                    fields::GW1z[i]=0.0;
                    fields::GW2x[i]=0.0;
                    fields::GW2y[i]=0.0;
                    fields::GW2z[i]=0.0;
                }
            }


            for(unsigned int i = 0 ; i < geom::ss ; i++)
            {
                const double s[3]={spins::sx(i),spins::sy(i),spins::sz(i)};
                double h[3]={fields::Hx(i),fields::Hy(i),fields::Hz(i)};

                //this is potentially slow because it is possible that the
                //compiler cannot optimize the calling of these functions
                //as they are pointers to functions defined at run time.
                fields::anisfp(h,s,i);
                fields::appliedfp(h);

                //exchange calculation
                for(unsigned int n = neigh::xadj[i] ; n < neigh::xadj[i+1] ; n++)
                {
                    const unsigned int neighbour = neigh::adjncy[n];
                    //look up the spin number
                    double sj[3]={spins::sx(neighbour),spins::sy(neighbour),spins::sz(neighbour)};
                    double area=neigh::surfArea[n];
                    double exchpf=2.0*tdp::sysexchstiff[i]/(tdp::sysme[i]*tdp::sysme[i]*mat::Ms);
                    //std::cout << exchpf/area << "\t" << area << "\t" << tdp::sysexchstiff[i] << "\t" << tdp::sysme[i] << "\t" << mat::Ms << std::endl;
                    //std::cin.get();
                    h[0]+=exchpf*((sj[0]-s[0])/area);
                    h[1]+=exchpf*((sj[1]-s[1])/area);
                    h[2]+=exchpf*((sj[2]-s[2])/area);

                }
                const double mi2=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
                const double oomi2=1./(mi2);

                //longitudinal field component
                if(tdp::systemp[i]<mat::Tc)
                {
                    const double oneover2chipar=0.5/(tdp::syschipar[i]);
                    const double bb=(1.0-mi2/(tdp::sysme[i]*tdp::sysme[i]));
                    h[0]+=oneover2chipar*bb*s[0];
                    h[1]+=oneover2chipar*bb*s[1];
                    h[2]+=oneover2chipar*bb*s[2];
                }
                else
                {
                    const double oneoverchipar=-1.0/tdp::syschipar[i];
                    const double bb=(1.0+3.0*mat::Tc*mi2/(5.0*(tdp::systemp[i]-mat::Tc)));
                    h[0]+=oneoverchipar*bb*s[0];
                    h[1]+=oneoverchipar*bb*s[1];
                    h[2]+=oneoverchipar*bb*s[2];
                }

                const double sdoth = s[0]*h[0]+s[1]*h[1]+s[2]*h[2];
                const double tpf=tdp::sysW2pf[i]*tdp::sysSigma[i];
                const double hpt[3]={h[0]+tpf*fields::GW2x[i],h[1]+tpf*fields::GW2y[i],h[2]+tpf*fields::GW2z[i]};
                const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
                const double sxhpt[3]={s[1]*hpt[2] - s[2]*hpt[1],s[2]*hpt[0]-s[0]*hpt[2],s[0]*hpt[1]-s[1]*hpt[0]};

    //            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};
                const double sxsxhpt[3]={s[1]*sxhpt[2]-s[2]*sxhpt[1],s[2]*sxhpt[0]-s[0]*sxhpt[2],s[0]*sxhpt[1]-s[1]*sxhpt[0]};
                const double aperp=tdp::sysalphaperp[i];
                const double apar=tdp::sysalphapar[i];
                const double tpf1=tdp::sysW1pf[i]*tdp::sysSigma[i];
                fnx[i] = sim::llbpf*(sxh[0] + (aperp*oomi2)*sxsxhpt[0] - (apar*oomi2)*sdoth*s[0]) + fields::GW1x[i]*tpf1;
                fny[i] = sim::llbpf*(sxh[1] + (aperp*oomi2)*sxsxhpt[1] - (apar*oomi2)*sdoth*s[1]) + fields::GW1y[i]*tpf1;
                fnz[i] = sim::llbpf*(sxh[2] + (aperp*oomi2)*sxsxhpt[2] - (apar*oomi2)*sdoth*s[2]) + fields::GW1z[i]*tpf1;

                spins::esx[i] = s[0] + fnx[i]*sim::rdt;
                spins::esy[i] = s[1] + fny[i]*sim::rdt;
                spins::esz[i] = s[2] + fnz[i]*sim::rdt;

            }

            for(unsigned int i = 0 ; i < geom::ss ; i++)
            {
                const double s[3]={spins::esx[i],spins::esy[i],spins::esz[i]};
                double h[3]={fields::Hx[i],fields::Hy[i],fields::Hz[i]};

                fields::anisfp(h,s,i);
                fields::appliedfp(h);

                //exchange calculation
                for(unsigned int n = neigh::xadj[i] ; n < neigh::xadj[i+1] ; n++)
                {
                    const unsigned int neighbour = neigh::adjncy[n];
                    //look up the spin number
                    double sj[3]={spins::esx(neighbour),spins::esy(neighbour),spins::esz(neighbour)};
                    double area=neigh::surfArea[n];
                    double exchpf=2.0*tdp::sysexchstiff[i]/(tdp::sysme[i]*tdp::sysme[i]*mat::Ms);
                    h[0]+=exchpf*((sj[0]-s[0])/area);
                    h[1]+=exchpf*((sj[1]-s[1])/area);
                    h[2]+=exchpf*((sj[2]-s[2])/area);

                }
                const double mi2=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
                const double oomi2=1./(mi2);

                //longitudinal field
                if(tdp::systemp[i]<mat::Tc)
                {
                    const double oneover2chipar=0.5/(tdp::syschipar[i]);
                    const double bb=(1.0-mi2/(tdp::sysme[i]*tdp::sysme[i]));
                    h[0]+=oneover2chipar*bb*s[0];
                    h[1]+=oneover2chipar*bb*s[1];
                    h[2]+=oneover2chipar*bb*s[2];
                }
                else
                {
                    const double oneoverchipar=-1.0/tdp::syschipar[i];
                    const double bb=(1.0+3.0*mat::Tc*mi2/(5.0*(tdp::systemp[i]-mat::Tc)));
                    h[0]+=oneoverchipar*bb*s[0];
                    h[1]+=oneoverchipar*bb*s[1];
                    h[2]+=oneoverchipar*bb*s[2];
                }

                const double tpf=tdp::sysW2pf[i]*tdp::sysSigma[i];
                const double hpt[3]={h[0]+tpf*fields::GW2x[i],h[1]+tpf*fields::GW2y[i],h[2]+tpf*fields::GW2z[i]};
                const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
                const double sxhpt[3]={s[1]*hpt[2] - s[2]*hpt[1],s[2]*hpt[0]-s[0]*hpt[2],s[0]*hpt[1]-s[1]*hpt[0]};
                const double sxsxhpt[3]={s[1]*sxhpt[2]-s[2]*sxhpt[1],s[2]*sxhpt[0]-s[0]*sxhpt[2],s[0]*sxhpt[1]-s[1]*sxhpt[0]};
    //            const double sxsxh[3]={s[1]*sxh[2]-s[2]*sxh[1],s[2]*sxh[0]-s[0]*sxh[2],s[0]*sxh[1]-s[1]*sxh[0]};
                const double sdoth = s[0]*h[0]+s[1]*h[1]+s[2]*h[2];
                const double aperp=tdp::sysalphaperp[i];
                const double apar=tdp::sysalphapar[i];
                const double tpf1=tdp::sysW1pf[i]*tdp::sysSigma[i];

                const double fnp1[3]={sim::llbpf*(sxh[0] + (aperp*oomi2)*sxsxhpt[0] - (apar*oomi2)*sdoth*s[0]) + fields::GW1x[i]*tpf1,
                                      sim::llbpf*(sxh[1] + (aperp*oomi2)*sxsxhpt[1] - (apar*oomi2)*sdoth*s[1]) + fields::GW1y[i]*tpf1,
                                      sim::llbpf*(sxh[2] + (aperp*oomi2)*sxsxhpt[2] - (apar*oomi2)*sdoth*s[2]) + fields::GW1z[i]*tpf1};

                spins::sx[i] = spins::sx[i] + 0.5*(fnx[i]+fnp1[0])*sim::rdt;
                spins::sy[i] = spins::sy[i] + 0.5*(fny[i]+fnp1[1])*sim::rdt;
                spins::sz[i] = spins::sz[i] + 0.5*(fnz[i]+fnp1[2])*sim::rdt;


            }
        }
        t+=fields::dfu;
    }
}
