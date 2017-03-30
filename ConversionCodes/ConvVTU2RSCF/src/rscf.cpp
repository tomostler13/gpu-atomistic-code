// File: discVTU.cpp
// Author:Tom Ostler
// Created: 25 Feb 2016
// Last-modified: 30 Mar 2017 10:34:28

//The purpose of this section of code is to take the output of the vtu
//files and put the spins into cells and calculate the average magnetization
//within the cell. This aids in visualising systems with too many spins.
//At the moment it only works for a single species because as of the creation
//date the VTU output format of the spin dynamics code does not have any
//knowledge of the species
#include <cmath>
#include <iostream>
#include <libconfig.h++>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <fftw3.h>
#include "../../../inc/error.h"
#include "../../../inc/array.h"
#include "../../../inc/array2d.h"
#include "../../../inc/array3d.h"
#include "../../../inc/array4d.h"
#include "../../../inc/defines.h"
#include "../../../inc/random.h"
#include "../../inc/inputs.h"
int main(int argc,char *argv[])
{
    Array4D<double> mag;
    //Array3D<int> coords;
    //Array2D<int> lu;
    Array2D<double> spins;
    Array3D<double> spinsx,spinsy,spinsz;
    Array3D<fftw_complex> kspinsx,kspinsy,kspinsz;
    fftw_plan SxP,SyP,SzP;
    fftw_plan SxB,SyB,SzB;
    int cplxdim=(inputs::dimin[2]/2)+1;
    std::cout << "#Reading config file..." << std::flush;
    inputs::readcff(argc,argv);

    double normfac=1./(static_cast<double>(inputs::dimin[0]*inputs::dimin[1]*inputs::dimin[2]));
    std::cout << "Done" << std::endl;

    std::cout << "#Resizing magnetization and spin arrays..." << std::flush;
    mag.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],3);
    mag.IFill(0);
    spins.resize(inputs::nspins,6);
    spinsx.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    spinsy.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    spinsz.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    spinsx.IFill(0);
    spinsy.IFill(0);
    spinsz.IFill(0);
    kspinsx.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    kspinsy.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    kspinsz.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    kspinsx.IFill(0);
    kspinsy.IFill(0);
    kspinsz.IFill(0);
    SxP=fftw_plan_dft_r2c_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],spinsx.ptr(),kspinsx.ptr(),FFTW_ESTIMATE);
    SyP=fftw_plan_dft_r2c_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],spinsy.ptr(),kspinsy.ptr(),FFTW_ESTIMATE);
    SzP=fftw_plan_dft_r2c_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],spinsz.ptr(),kspinsz.ptr(),FFTW_ESTIMATE);
    SxB=fftw_plan_dft_c2r_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],kspinsx.ptr(),spinsx.ptr(),FFTW_ESTIMATE);
    SyB=fftw_plan_dft_c2r_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],kspinsy.ptr(),spinsy.ptr(),FFTW_ESTIMATE);
    SzB=fftw_plan_dft_c2r_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],kspinsz.ptr(),spinsz.ptr(),FFTW_ESTIMATE);
    std::cout << "Done" << std::endl;

    std::cout << "#inputs::dimin = " << inputs::dimin[0] << " , " << inputs::dimin[1] << " , " << inputs::dimin[2] << std::endl;
    for(unsigned int t = inputs::lt ; t < inputs::ut+1 ; t+=inputs::dt)
    {
        std::stringstream sstr;
        sstr << inputs::fb << t << ".vtu";
        std::string str=sstr.str();

        std::cout << "#Opening file " << str << "..." << std::flush;
        std::ifstream ifs(str.c_str());
        if(!ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open input file");
        }
        std::cout << "done" << std::endl;
        std::string dumpline;
        //get rid of the top six lines
        for(unsigned int i = 0 ; i < 6 ; i++)
        {
            std::getline(ifs,dumpline);
//            std::cout << dumpline << std::endl;
        }
        //get the spin vectors (sx,sy,sz)
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            ifs >> spins(i,0) >> spins(i,1) >> spins(i,2);
//            std::cout << spins(i,0) << "\t" << spins(i,1) << "\t" << spins(i,2) << std::endl;
        }
//        std::cout << spins(inputs::nspins-1,0) << std::endl;
        //get rid of the middle six lines
        for(unsigned int i = 0 ; i < 7 ; i++)
        {
            std::getline(ifs,dumpline);
//            std::cout << dumpline << std::endl;
        }
//        std::cin.get();
        //get the coordinates
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            ifs >> spins(i,3) >> spins(i,4) >> spins(i,5);
            int ci[3]={int((spins(i,3)/inputs::abc[0])+0.5),int((spins(i,4)/inputs::abc[1])+0.5),int((spins(i,5)/inputs::abc[2])+0.5)};
            spinsx(ci[0],ci[1],ci[2])=spins(i,0);
            spinsy(ci[0],ci[1],ci[2])=spins(i,1);
            spinsz(ci[0],ci[1],ci[2])=spins(i,2);
            //std::cout << ci[0] << "\t" << ci[1] << "\t" << ci[2] << std::endl;
            //std::cout << spins(i,0) << "\t" << spins(i,1) << "\t" << spins(i,2) << std::endl;
        }
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
        for(unsigned int i = 0 ; i < inputs::dimin[0] ; i++)
        {
            for(unsigned int j = 0 ; j < inputs::dimin[1] ; j++)
            {
                for(unsigned int k = 0 ; k < cplxdim ; k++)
                {
                    double tx[2]={0,0};
                    tx[0]=kspinsx(i,j,k)[0]*kspinsx(i,j,k)[0]+kspinsx(i,j,k)[1]*kspinsx(i,j,k)[1];
                    tx[1]=kspinsx(i,j,k)[1]*kspinsx(i,j,k)[0]-kspinsx(i,j,k)[0]*kspinsx(i,j,k)[1];
                    kspinsx(i,j,k)[0]=tx[0];
                    kspinsx(i,j,k)[1]=tx[1];
                    double ty[2]={0,0};
                    ty[0]=kspinsy(i,j,k)[0]*kspinsy(i,j,k)[0]+kspinsy(i,j,k)[1]*kspinsy(i,j,k)[1];
                    ty[1]=kspinsy(i,j,k)[1]*kspinsy(i,j,k)[0]-kspinsy(i,j,k)[0]*kspinsy(i,j,k)[1];
                    kspinsy(i,j,k)[0]=ty[0];
                    kspinsy(i,j,k)[1]=ty[1];
                    double tz[2]={0,0};
                    tz[0]=kspinsz(i,j,k)[0]*kspinsz(i,j,k)[0]+kspinsz(i,j,k)[1]*kspinsz(i,j,k)[1];
                    tz[1]=kspinsz(i,j,k)[1]*kspinsz(i,j,k)[0]-kspinsz(i,j,k)[0]*kspinsz(i,j,k)[1];
                    kspinsz(i,j,k)[0]=tz[0];
                    kspinsz(i,j,k)[1]=tz[1];
                }
            }
        }
        fftw_execute(SxB);
        fftw_execute(SyB);
        fftw_execute(SzB);
        ifs.close();
        if(ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not close the ifs file");
        }
        std::stringstream osstr;
        osstr << "rscf_" << t << ".dat";
        std::string ostr=osstr.str();
        std::ofstream ofs;
        ofs.open(ostr.c_str());
        if(!ofs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Couldn't open an output file");
        }
        //output in 2D (x,y hard coded)
        for(int i = 0 ; i < inputs::dimin[0] ; i++)
        {
            int ri=i-(inputs::dimin[0]/2);
            int li=0;
            if(ri<0)
            {
                li=ri+inputs::dimin[0];
            }
            else
            {
                li=ri;
            }
            for(int j = 0 ; j < inputs::dimin[1] ; j++)
            {
                int rj=j-(inputs::dimin[1]/2);
                int lj=0;
                if(rj<0)
                {
                    lj=rj+inputs::dimin[1];
                }
                else
                {
                    lj=rj;
                }
                //output sxsx, sysy and szsz correlation function
                ofs << ri << "\t" << rj << "\t" << spinsx(li,lj,0)*normfac << "\t" << spinsy(li,lj,0)*normfac << "\t" << spinsz(li,lj,0)*normfac << std::endl;
            }
        }
        ofs.close();
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
    }

    return(EXIT_SUCCESS);
}
