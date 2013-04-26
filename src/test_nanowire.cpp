// File: test_nanowire.cpp
// Author: Tom Ostler
// Created: 10 April 2013
// Last-modified: 25 Apr 2013 19:51:45
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>
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
    //fields of cold, hot and wire
    double CF=0.0,HF=0.0,WF=0.0;
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
    setting.lookupValue("Hot_Field",HF);
    setting.lookupValue("Cold_Field",CF);
    setting.lookupValue("Wire_Field",WF);
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
    FIXOUT(config::Info,"Field in cold part:" << CF << std::endl);
    FIXOUT(config::Info,"Field in hot part:" << HF << std::endl);
    FIXOUT(config::Info,"Field in wire:" << WF << std::endl);
    FIXOUT(config::Info,"Equilibration time:" << et << std::endl);
    FIXOUT(config::Info,"Run time:" << rt << std::endl);
    FIXOUT(config::Info,"Number of equilibration timesteps:" << ets << std::endl);
    FIXOUT(config::Info,"Number of runtime steps:" << rts << std::endl);
    FIXOUT(config::Info,"Setting temperature:" << std::flush);
    //counters for number of hot, cold and wire spins
    int cc=0,wc=0,hc=0,maxz;
    for(unsigned int i = 0 ; i < geom::nspins ; i++)
    {
        fields::HAppy[i]=0;
        fields::HAppz[i]=0;

        int coords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
        if(coords[2]>maxz)
        {
            maxz=coords[2];
        }
        if(coords[2]<geom::cut0)
        {
            llg::osT[i]=LTR;
            fields::HAppx[i]=CF;
            cc++;
        }
        else if(coords[2]>=geom::cut0 && coords[2]<geom::cut1)// < geom::cut1 && coords[0] > (geom::dim[0]-geom::width)/2 && coords[0] < (geom::dim[0]+geom::width)/2 && coords[1] > (geom::dim[1]-geom::width)/2 && coords[1] < (geom::dim[1]+geom::width)/2)
        {
            llg::osT[i]=TOW;
            fields::HAppx[i]=WF;
            wc++;
        }
        else if(coords[2]>=geom::cut1)
        {
            hc++;
            llg::osT[i]=HTR;
            fields::HAppx[i]=HF;
        }
    }
    //at either end of the wire we want to set the damping to critical to try to dampen
    //spinwaves so they are not reflected
    unsigned int la=geom::coords(0,0,0,0),ra=geom::coords(0,0,maxz,0);

    mat::lambda[la]=1.0;
    mat::lambda[ra]=1.0;
    mat::sigma[la]=sqrt(2.0*1.38e-23*mat::lambda[la]/(mat::mu[la]*mat::muB*llg::dt*mat::gamma));
    mat::sigma[ra]=sqrt(2.0*1.38e-23*mat::lambda[ra]/(mat::mu[ra]*mat::muB*llg::dt*mat::gamma));
    llg::llgpf[la] = -1./(1.0+mat::lambda[la]*mat::lambda[la]);
    llg::llgpf[ra] = -1./(1.0+mat::lambda[ra]*mat::lambda[ra]);
    std::ofstream ws,cs,hs;
    cs.open("coldSSF.dat");
    ws.open("wireSSF.dat");
    hs.open("hotSSF.dat");
    if(!cs.is_open() || !ws.is_open() || ! hs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file for outputting structure factors");
    }
    SUCCESS(config::Info);
    Array<fftw_complex> SpmCold(cc),SpmHot(hc),SpmWire(wc);
    fftw_plan cp,hp,wp;
    cp=fftw_plan_dft_1d(cc,SpmCold.ptr(),SpmCold.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
    hp=fftw_plan_dft_1d(hc,SpmHot.ptr(),SpmHot.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
    wp=fftw_plan_dft_1d(wc,SpmWire.ptr(),SpmWire.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
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
    std::ofstream mx_chain_out("mx_time.dat");
    mx_chain_out << "#\t" << std::flush;
    int points_of_interest[9]={1,geom::cut0/2,geom::cut0-1,geom::cut0,(geom::cut1-geom::cut0)/2+geom::cut0,geom::cut1-1,geom::cut1,geom::cut1+(geom::dim[2]-geom::cut1)/2,geom::dim[2]-1};
    unsigned int count=0;
    for(unsigned int i = 0 ; i < geom::nspins ; i++)
    {
        int coords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
        if(coords[2]==points_of_interest[count])
        {
            mx_chain_out << points_of_interest[count] << "\t" << std::flush;
            count++;
        }
    }
    mx_chain_out << std::endl;
    bool topcheck=false;
    for(unsigned int t = ets ; t < rts+ets ; t++)
    {
        llg::integrate(t);
        if(t%spins::update==0)
        {
            count=0;
            unsigned int ccount=0,hcount=0,wcount=0;
            mx_chain_out << t << "\t";
            for(int i = 0 ; i < geom::nspins ; i++)
            {
                int coords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                if(coords[2]==points_of_interest[count])
                {
                    mx_chain_out << spins::Sx[i] << "\t";
                    count++;
                }
                if(coords[2]<geom::cut0)
                {
                    SpmCold(ccount)[0]=spins::Sx[i];
                    SpmCold(ccount)[1]=spins::Sy[i];
                    ccount++;
                }
                else if(coords[2]>=geom::cut0 && coords[2]<geom::cut1)
                {
                    SpmWire(wcount)[0]=spins::Sx[i];
                    SpmWire(wcount)[1]=spins::Sy[i];
                    wcount++;
                }
                else if(coords[2]>=geom::cut1)
                {
                    SpmHot(hcount)[0]=spins::Sx[i];
                    SpmHot(hcount)[1]=spins::Sy[i];
                    hcount++;
                }
            }
            mx_chain_out << std::endl;
            if(ccount!=cc)
            {
                std::cerr << "BLA1" << std::endl;
            }
            if(hcount!=hc)
            {
                std::cerr << "BLA2" << std::endl;
            }
            if(wcount!=wc)
            {
                std::cerr << "BLA3" << std::endl;
            }

            fftw_execute(cp);
            fftw_execute(hp);
            fftw_execute(wp);
            //output wire structure factor
            for(int i = wc/2 ; i < wc ; i++)
            {
                ws << t << "\t" << i-wc << "\t" << sqrt(SpmWire(i)[0]*SpmWire(i)[0]+SpmWire(i)[1]*SpmWire(i)[1]) << std::endl;
            }
            for(int i = 1 ; i < wc/2 ; i++)
            {
                ws << t << "\t" << i << "\t" << sqrt(SpmWire(i)[0]*SpmWire(i)[0]+SpmWire(i)[1]*SpmWire(i)[1]) << std::endl;
            }
            ws << std::endl;
            //output cold structure factor
            for(int i = cc/2 ; i < cc ; i++)
            {
                cs <<  t << "\t" << i-cc << "\t" << sqrt(SpmCold(i)[0]*SpmCold(i)[0]+SpmCold(i)[1]*SpmCold(i)[1]) << std::endl;
            }
            for(int i = 1 ; i < cc/2 ; i++)
            {
                cs << t << "\t" << i << "\t" << sqrt(SpmCold(i)[0]*SpmCold(i)[0]+SpmCold(i)[1]*SpmCold(i)[1]) << std::endl;
            }
            cs << std::endl;
            //output hot structure factor
            for(int i = hc/2 ; i < hc ; i++)
            {
                hs << t << "\t" << i-hc << "\t" << sqrt(SpmHot(i)[0]*SpmHot(i)[0]+SpmHot(i)[1]*SpmHot(i)[1]) << std::endl;
            }
            for(int i = 1 ; i < cc/2 ; i++)
            {
                hs << t << "\t" << i << "\t" << sqrt(SpmHot(i)[0]*SpmHot(i)[0]+SpmHot(i)[1]*SpmHot(i)[1]) << std::endl;
            }
            hs << std::endl;
            const double mx = util::reduceCPU(spins::Sx,geom::nspins)/double(geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins)/double(geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins)/double(geom::nspins);
            const double modm=sqrt(mx*mx+my*my+mz*mz);
            magout << t << "\t" << mx << "\t" << my << "\t" << mz << "\t" << modm << std::endl;
            //util::outputSpinsVTU(t);
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
    mx_chain_out.close();



}
