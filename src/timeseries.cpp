// File: timeseries.h
// Author: Tom Ostler
// Created: 05 April 2013
// Last-modified: 07 Apr 2013 21:23:55
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
void sim::timeseries(int argc,char *argv[])
{
    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Timeseries details***" << std::endl;
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

    libconfig::Setting &setting = config::cfg.lookup("timeseries");
    double et,rt,temperature=0.0;
    unsigned int ets,rts;
    setting.lookupValue("EquilTime",et);
    setting.lookupValue("RunTime",rt);
    FIXOUT(config::Info,"Equilibration time: " << et << std::endl);
    FIXOUT(config::Info,"Run time:" << rt << std::endl);
    ets = int(et/llg::dt);
    rts = int(rt/llg::dt);
    FIXOUT(config::Info,"Equilibration timesteps:" << ets << std::endl);
    FIXOUT(config::Info,"Run timesteps:" << rts << std::endl);
    setting.lookupValue("Temperature",temperature);
    FIXOUT(config::Info,"Temperature:" << temperature << std::endl);

    bool outmag=false;
    std::ofstream magout,emagout;
    std::string MagFilestr;
    setting.lookupValue("outmag",outmag);
    FIXOUT(config::Info,"Output magnetization data for time series to:" << config::isTF(outmag) << std::endl);
    if(outmag==true)
    {
        setting.lookupValue("MagFileName",MagFilestr);
        FIXOUT(config::Info,"Magnetization files output to:" << MagFilestr << std::endl);
    }
    if(outmag)
    {
        std::stringstream MagFilesstr,eqsstr;
        MagFilesstr << MagFilestr << ".dat";
        eqsstr << MagFilestr << "_eq" << ".dat";
        std::string mopf=MagFilesstr.str(),eopf=eqsstr.str();
        magout.open(mopf.c_str());
        emagout.open(eopf.c_str());
        if(!magout.is_open() || !emagout.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for outputting magnetization data");
        }
    }
    llg::T=temperature;
    //work out size for 4D time series data
    unsigned int nt=rts/spins::update;
    double t_sample=double(spins::update)*llg::dt;
    FIXOUT(config::Info,"Total number of time samples:" << nt << " timesteps" << std::endl);
    FIXOUT(config::Info,"Temporal sampling (real time)" << spins::update*llg::dt << " seconds" << std::endl);
//    Array4D<fftw_complex> tsd;
//    tsd.IFill(0);
//    tsd.resize(nt,geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
    unsigned int tempcounter=0;
    for(unsigned int t = 0 ; t < ets ; t++)
    {
        llg::integrate(t);
        if(t%spins::update==0 && outmag)
        {
            const double mx = util::reduceCPU(spins::Sx,geom::nspins);
            const double my = util::reduceCPU(spins::Sy,geom::nspins);
            const double mz = util::reduceCPU(spins::Sz,geom::nspins);
            emagout << double(t)*llg::dt << "\t" << mx/double(geom::nspins) << "\t" << my/double(geom::nspins) << "\t" << mz/double(geom::nspins) << std::endl;
        }
    }
    emagout.close();
    if(emagout.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close equilibration magnetisation file");
    }
    std::ofstream binfile("spins_ts.dat",std::ofstream::binary);
    if(binfile.is_open()!=true)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file for writing binary information");
    }
    fftw_complex *fftw_data=NULL;
    try
    {
        fftw_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*geom::nspins);
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
            if(outmag)
            {
                const double mx = util::reduceCPU(spins::Sx,geom::nspins);
                const double my = util::reduceCPU(spins::Sy,geom::nspins);
                const double mz = util::reduceCPU(spins::Sz,geom::nspins);
                magout << double(t)*llg::dt << "\t" << mx/double(geom::nspins) << "\t" << my/double(geom::nspins) << "\t" << mz/double(geom::nspins) << std::endl;
            }
            /*for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                unsigned int coords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                if(tempcounter<nt)
                {
                    //tsd(tempcounter,coords[0],coords[1],coords[2])[0]=spins::Sx(i);
                    //tsd(tempcounter,coords[0],coords[1],coords[2])[1]=spins::Sy(i);
                    tempcounter++;
                }
            }*/

            unsigned int atom_counter=0;
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
                            fftw_data[atom_counter][0]=spins::Sx(sn);
                            fftw_data[atom_counter][1]=spins::Sy(sn);
//                            std::cout << fftw_data[atom_counter][0] << "\t" << fftw_data[atom_counter][1] << "\t" << sn << std::endl;
                            atom_counter++;
                        }
                    }
                }
            }
            binfile.write(reinterpret_cast<char*>(fftw_data),geom::nspins*sizeof(fftw_complex));
        }
    }
    try
    {
        fftw_free(fftw_data);
    }
    catch(...)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not free fftw data");
    }
    binfile.close();
    if(binfile.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close binfile");
    }
    //int n[4]={nt,geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]};
    //int nk[4]={n[0]/2,n[1]/2,n[2]/2,n[3]/2};
//    double sotsd=double(n[0]*n[1]*n[2]*n[3]);
//    FIXOUT(config::Info,"Planning 4D fourier transform:" << std::flush);
//    fftw_plan plan4d = fftw_plan_dft(4,n,tsd.ptr(),tsd.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
//    SUCCESS(config::Info);
//    FIXOUT(config::Info,"Performing 4D fourier transform:" << std::flush);
//    fftw_execute(plan4d);
/*    SUCCESS(config::Info);

    std::stringstream filess;
 std::cout<<"Printing files..."<<std::endl;
    for(unsigned int i=0;i<n[1];++i){
      for(unsigned int j=0;j<n[2];++j){
        for(unsigned int k=0;k<n[3];++k){
        if(i==0 && j==0)
        {
          filess<<"kx"<<i<<"ky"<<j<<"kz"<<k<<".dat";
          std::string filename = filess.str();
          std::ofstream outfile (filename.c_str());
          for(int t=0;t<nk[0];++t){
          const double freq = (t*2.0*M_PI)/(t_sample*double(n[0]));
          // 4D array lookup (sic!)
          unsigned int idx1 = k+n[3]*(j+n[2]*(i+n[1]*t));
          unsigned int idx2 = k+n[3]*(j+n[2]*(i+n[1]*(n[0]-t-1)));
          fftw_complex p1 = {tsd(t,i,j,k)[0],tsd(t,i,j,k)[1]};
          fftw_complex p2 = {tsd(n[0]-t-1,i,j,k)[0],tsd(n[0]-t-1,i,j,k)[1]};
          const double a = p1[0]*p1[0]+p1[1]*p1[1];
          const double b = p2[0]*p2[0]+p2[1]*p2[1];
          outfile <<t<<"\t"<< freq <<"\t"<<a<<"\n";
          outfile <<-t<<"\t"<< -freq <<"\t"<<b<<"\n";
        }
        outfile.close();
        filess.str("");
        }
      }
    }
  }
  std::cout<<"Done"<<std::endl;
  */
/*    for(int i = -nk[1] ; i < nk[1]-1 ; i++)
    {
        for(int j = -nk[2] ; j < nk[2]-1 ; j++)
        {
            for(int k = -nk[3] ; k < nk[3]-1 ; k++)
            {
                if(i==0 && j==0)
                {
                //for each q-vector output the frequency
                std::stringstream sstr;
                sstr << "kx_" << i << "_ky_" << j << "_kz_" << k << ".dat";
                std::string str=sstr.str();
                std::ofstream opf(str.c_str());
                std::streambuf *buf=std::cout.rdbuf(opf.rdbuf());
                int luc[3]={i,j,k};
                for(unsigned int check=0 ; check < 3 ; check++)
                {
                    if(luc[check]<0)
                    {
                        luc[check]=geom::dim[0]*geom::Nk[0]+luc[check];
                    }
                }
                for(int t=-nk[0] ; t < nk[0]-1 ; t++)
                {
                    int lut=t;
                    if(lut < 0)
                    {
                        lut=nt+lut;
                    }
                    double freq=double(t)*2.0*M_PI/(rt);
                    fftw_complex p1={tsd(lut,luc[0],luc[1],luc[2])[0],tsd(lut,luc[0],luc[1],luc[2])[1]};
                    std::cout << t << "\t" << freq << "\t" << p1[0]*p1[0]+p1[1]*p1[1] << std::endl;
                }
                opf.close();
                }
            }
        }
    }*/
}
