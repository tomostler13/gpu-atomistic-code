// File: test_nanowire.cpp
// Author: Tom Ostler
// Created: 10 April 2013
// Last-modified: 19 Apr 2013 13:00:48
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "../inc/arrays.h"
#include "../inc/defines.h"
#include "../inc/util.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/intmat.h"
#include "../inc/fields.h"
#include "../inc/llg.h"
#include "../inc/sim.h"
#include "../inc/exch.h"
void sim::test_nanowire(int argc,char *argv[])
{
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Test Nanowire details***" << std::endl;
    try
    {
        config::cfg.readFile(argv[1]);
    }
    catch(const libconfig::FileIOException &fioex)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("I/O error while reading config file");
    }
    catch(const libconfig::ParseException &pex)
    {
        error::errPreamble(__FILE__,__LINE__);
        std::cerr << ". Parse error at " << pex.getFile()  << ":" << pex.getLine() << "-" << pex.getError() << "***\n" << std::endl;
        exit(EXIT_FAILURE);
    }

    libconfig::Setting &setting = config::cfg.lookup("test_nanowire");

    //low temperature reservoir
    double LTR=0.0;
    //high temperature reservoir
    double HTR=0.0;
    //temperature of wire
    double TOW=0.0;
    //equilibration and runtime
    double et=0.0,rt=0.0;
    //steps
    bool outmag=false;
    std::string MagFilestr;
    setting.lookupValue("EquilTime",et);
    setting.lookupValue("RunTime",rt);
    setting.lookupValue("Wire_Temperature",TOW);
    setting.lookupValue("High_Temperature",HTR);
    setting.lookupValue("Low_Temperature",LTR);
    setting.lookupValue("outmag",outmag);

    unsigned int ets=int(et/llg::dt),rts=int(rt/llg::dt);
    FIXOUT(config::Info,"Output magnetization data for each temperature:" << config::isTF(outmag) << std::endl);
    if(outmag==true)
    {
        setting.lookupValue("MagFileName",MagFilestr);
        FIXOUT(config::Info,"Magnetization files output to:" << MagFilestr << std::endl);
    }
    std::ofstream magout;
    if(outmag)
    {
        std::stringstream MagFilesstr,eqsstr;
        MagFilesstr << MagFilestr << "_" << LTR << "_" << TOW << "_" << HTR << ".dat";
        std::string mopf=MagFilesstr.str();
        magout.open(mopf.c_str());
        if(!magout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for outputting magnetization data");
        }
    }

    FIXOUT(config::Info,"Temperature of high temperature reservoir:" << HTR << std::endl);
    FIXOUT(config::Info,"Temperature of low temperature reservoir:" << LTR << std::endl);
    FIXOUT(config::Info,"Temperature of wire:" << TOW << std::endl);
    FIXOUT(config::Info,"Equilibration time:" << et << std::endl);
    FIXOUT(config::Info,"Run time:" << rt << std::endl);
    FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    FIXOUT(config::Info,"Number of runtime steps:" << rts << std::endl);
    FIXOUT(config::Info,"Setting temperature:" << std::flush);
    //counters for number of hot, cold and wire spins
    unsigned int cc=0,wc=0,hc=0;
    for(unsigned int i = 0 ; i < geom::nspins ; i++)
    {
        int coords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
        if(coords[2]<geom::cut0)
        {
            llg::osT[i]=LTR;
            cc++;
//            std::cout << "BALSDJKLASDFJKLSADF" << std::endl;
        }
        else if(coords[2]>=geom::cut0 && coords[2]<geom::cut1)// < geom::cut1 && coords[0] > (geom::dim[0]-geom::width)/2 && coords[0] < (geom::dim[0]+geom::width)/2 && coords[1] > (geom::dim[1]-geom::width)/2 && coords[1] < (geom::dim[1]+geom::width)/2)
        {
            llg::osT[i]=TOW;
            wc++;
//            std::cout << "________________________________" << std::endl;
        }
        else if(coords[2]>=geom::cut1)
        {
            hc++;
            llg::osT[i]=HTR;
//            std::cout << "Setting" << std::endl;
        }
    }

    SUCCESS(config::Info);
    FIXOUT(config::Info,"Number of atoms in cold part:" << cc << std::endl);
    FIXOUT(config::Info,"Number of atoms in wire:" << wc << std::endl);
    FIXOUT(config::Info,"Number of atoms in hot part:" << hc << std::endl);
    if(llg::rscf)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("You have set output of correlation function to true when you do not have translational invariance, setting to false.");
        llg::rscf=false;
    }
    for(unsigned int t = 0 ; t < ets ; t++)
    {
        llg::integrate(t);
        if(t%spins::update==0)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);
            magout << "\t" << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            //util::outputSpinsVTU(t);
        }
    }
    std::ofstream binfilew("spins_ts_wire.dat",std::ofstream::binary);
    std::ofstream binfilec("spins_ts_cold.dat",std::ofstream::binary);
    std::ofstream binfileh("spins_ts_hot.dat",std::ofstream::binary);
    if(binfilew.is_open()!=true || binfilec.is_open()!=true || binfileh.is_open()!=true)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file for writing binary information");
    }
    fftw_complex *fftw_dataw=NULL;
    fftw_complex *fftw_datac=NULL;
    fftw_complex *fftw_datah=NULL;
    try
    {
        fftw_dataw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*wc);
        fftw_datac = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*cc);
        fftw_datah = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*hc);
    }
    catch(...)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not malloc fftw array");
    }
    for(unsigned int t = ets ; t < rts+ets ; t++)
    {
        llg::integrate(t);
        if(t%spins::update==0)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);
            magout << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            util::outputSpinsVTU(t);
            unsigned int wire_counter=0,cold_counter=0,hot_counter=0;
            for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                    {

                        int sn=geom::coords(i,j,k,0);
    //                    std::cout << sn << std::endl;
                        if(sn>-1)
                        {
                            if(k<geom::Nk[2]*geom::cut0)
                            {
                                fftw_datac[cold_counter][0]=spins::Sx(sn);
                                fftw_datac[cold_counter][1]=spins::Sy(sn);
                                cold_counter++;
                            }
                            else if(k>=geom::Nk[2]*geom::cut0 && k < geom::Nk[2]*geom::cut1)
                            {
                                fftw_dataw[wire_counter][0]=spins::Sx(sn);
                                fftw_dataw[wire_counter][1]=spins::Sy(sn);
                                wire_counter++;
                            }
                            else if(k>=geom::Nk[2]*geom::cut1)
                            {
                                fftw_datah[hot_counter][0]=spins::Sx(sn);
                                fftw_datah[hot_counter][1]=spins::Sy(sn);
                                hot_counter++;
                            }
                        }
                    }
                }
            }
            if(wire_counter!=wc || cold_counter!=cc || hot_counter!=hc)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Incorrect number of spins found in wire");
            }
            binfilew.write(reinterpret_cast<char*>(fftw_dataw),wc*sizeof(fftw_complex));
            binfilec.write(reinterpret_cast<char*>(fftw_datac),cc*sizeof(fftw_complex));
            binfileh.write(reinterpret_cast<char*>(fftw_datah),hc*sizeof(fftw_complex));

        }
    }



}
