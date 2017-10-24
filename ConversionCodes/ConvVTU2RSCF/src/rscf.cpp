// File: discVTU.cpp
// Author:Tom Ostler
// Created: 25 Feb 2016
// Last-modified: 31 Mar 2017 15:33:11

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
    Array3D<double> spinsx,spinsy,spinsz,spinsmod,Cxx,Cyy,Czz,Csds;
    Array3D<fftw_complex> kspinsx,kspinsy,kspinsz,kspinsmod,CKxx,CKyy,CKzz,CKsds;
    fftw_plan SxP,SyP,SzP;
    fftw_plan SxB,SyB,SzB,SdotB;
    inputs::readcff(argc,argv);
    inputs::dimin[2]+=1;
    std::cout << "#Reading config file..." << std::flush;

    int cplxdim=(inputs::dimin[2]/2)+1;


    double normfac=(1./(static_cast<double>(inputs::dimin[0]*inputs::dimin[1]*inputs::dimin[2])))*(1./(static_cast<double>(inputs::dimin[0]*inputs::dimin[1]*inputs::dimin[2])));
    std::cout << "Done" << std::endl;

    std::cout << "New dim size " << inputs::dimin[2] << std::endl;

    std::cout << "#Normalisation factor" << normfac << std::endl;

    std::cout << "#Resizing magnetization and spin arrays..." << std::flush;
    mag.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],3);
    mag.IFill(0);
    spins.resize(inputs::nspins,6);
    spinsx.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    spinsy.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    spinsz.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    spinsmod.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    Cxx.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    Cyy.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    Czz.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    Csds.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);
    CKxx.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    CKyy.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    CKzz.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    CKsds.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    spinsx.IFill(0);
    spinsy.IFill(0);
    spinsz.IFill(0);
    spinsmod.IFill(0);
    Cxx.IFill(0);
    Cyy.IFill(0);
    Czz.IFill(0);
    Csds.IFill(0);
    CKxx.IFill(0);
    CKyy.IFill(0);
    CKzz.IFill(0);
    CKsds.IFill(0);
    kspinsx.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    kspinsy.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    kspinsz.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    kspinsmod.resize(inputs::dimin[0],inputs::dimin[1],cplxdim);
    kspinsx.IFill(0);
    kspinsy.IFill(0);
    kspinsz.IFill(0);
    kspinsmod.IFill(0);
    SxP=fftw_plan_dft_r2c_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],spinsx.ptr(),kspinsx.ptr(),FFTW_ESTIMATE);
    SyP=fftw_plan_dft_r2c_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],spinsy.ptr(),kspinsy.ptr(),FFTW_ESTIMATE);
    SzP=fftw_plan_dft_r2c_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],spinsz.ptr(),kspinsz.ptr(),FFTW_ESTIMATE);
    SxB=fftw_plan_dft_c2r_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],CKxx.ptr(),Cxx.ptr(),FFTW_ESTIMATE);
    SyB=fftw_plan_dft_c2r_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],CKyy.ptr(),Cyy.ptr(),FFTW_ESTIMATE);
    SzB=fftw_plan_dft_c2r_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],CKzz.ptr(),Czz.ptr(),FFTW_ESTIMATE);
    SdotB=fftw_plan_dft_c2r_3d(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],CKsds.ptr(),Csds.ptr(),FFTW_ESTIMATE);
    spinsx.IFill(0);
    spinsy.IFill(0);
    spinsz.IFill(0);
    spinsmod.IFill(0);
    kspinsx.IFill(0);
    kspinsy.IFill(0);
    kspinsz.IFill(0);
    kspinsmod.IFill(0);
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
        spins.IFill(0);
        spinsx.IFill(0);
        spinsy.IFill(0);
        spinsz.IFill(0);
        kspinsx.IFill(0);
        kspinsy.IFill(0);
        kspinsz.IFill(0);
        kspinsmod.IFill(0);
        spinsmod.IFill(0);
        //get the spin vectors (sx,sy,sz)
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            ifs >> spins(i,0) >> spins(i,1) >> spins(i,2);
//            const double n=sqrt(spins(i,0)*spins(i,0)+spins(i,1)*spins(i,1)+spins(i,2)*spins(i,2));
//            spins(i,0)/=n;
//            spins(i,1)/=n;
//            spins(i,2)/=n;
            //std::cout << spins(i,0) << "\t" << spins(i,1) << "\t" << spins(i,2) << std::endl;
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
                    CKxx(i,j,k)[0]=kspinsx(i,j,k)[0]*kspinsx(i,j,k)[0]+kspinsx(i,j,k)[1]*kspinsx(i,j,k)[1];
                    CKxx(i,j,k)[1]=kspinsx(i,j,k)[1]*kspinsx(i,j,k)[0]-kspinsx(i,j,k)[0]*kspinsx(i,j,k)[1];
                    CKyy(i,j,k)[0]=kspinsy(i,j,k)[0]*kspinsy(i,j,k)[0]+kspinsy(i,j,k)[1]*kspinsy(i,j,k)[1];
                    CKyy(i,j,k)[1]=kspinsy(i,j,k)[1]*kspinsy(i,j,k)[0]-kspinsy(i,j,k)[0]*kspinsy(i,j,k)[1];
                    CKzz(i,j,k)[0]=kspinsz(i,j,k)[0]*kspinsz(i,j,k)[0]+kspinsz(i,j,k)[1]*kspinsz(i,j,k)[1];
                    CKzz(i,j,k)[1]=kspinsz(i,j,k)[1]*kspinsz(i,j,k)[0]-kspinsz(i,j,k)[0]*kspinsz(i,j,k)[1];

                    CKsds(i,j,k)[0]=kspinsx(i,j,k)[0]*kspinsx(i,j,k)[0] + kspinsx(i,j,k)[1]*kspinsx(i,j,k)[1] + kspinsy(i,j,k)[0]*kspinsy(i,j,k)[0]+kspinsy(i,j,k)[1]*kspinsy(i,j,k)[1] + kspinsz(i,j,k)[0]*kspinsz(i,j,k)[0]+kspinsz(i,j,k)[1]*kspinsz(i,j,k)[1];
                    CKsds(i,j,k)[1]=kspinsx(i,j,k)[1]*kspinsx(i,j,k)[0] - kspinsx(i,j,k)[0]*kspinsx(i,j,k)[1] + kspinsy(i,j,k)[1]*kspinsy(i,j,k)[0]-kspinsy(i,j,k)[0]*kspinsy(i,j,k)[1] + kspinsz(i,j,k)[1]*kspinsz(i,j,k)[0]-kspinsz(i,j,k)[0]*kspinsz(i,j,k)[1];
                }
            }
        }
        fftw_execute(SxB);
        fftw_execute(SyB);
        fftw_execute(SzB);
        fftw_execute(SdotB);
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
                ofs << ri << "\t" << rj << "\t" << Cxx(li,lj,0)*normfac << "\t" << Cyy(li,lj,0)*normfac << "\t" << Czz(li,lj,0)*normfac << "\t" << Csds(li,lj,0)*normfac << "\t" << sqrt(CKxx(li,lj,0)[0]*CKxx(li,lj,0)[0]+CKxx(li,lj,0)[1]*CKxx(li,lj,0)[1])*normfac << "\t" << sqrt(CKyy(li,lj,0)[0]*CKyy(li,lj,0)[0]+CKyy(li,lj,0)[1]*CKyy(li,lj,0)[1])*normfac << "\t" << sqrt(CKzz(li,lj,0)[0]*CKzz(li,lj,0)[0]+CKzz(li,lj,0)[1]*CKzz(li,lj,0)[1])*normfac << "\t" << sqrt(CKsds(li,lj,0)[0]*CKsds(li,lj,0)[0]+CKsds(li,lj,0)[1]*CKsds(li,lj,0)[1])*normfac << std::endl;
            }
            ofs << std::endl << std::endl;
        }
        /*for(int i = inputs::dimin[0]/2 ; i < inputs::dimin[0] ; i++)
        {
            ofs << i-static_cast<int>(inputs::dimin[0]) << "\t" << Csds(i,0,0)*normfac << std::endl;
        }
        for(int i = 0 ; i < inputs::dimin[0]/2 ; i++)
        {
            ofs << i << "\t" << Csds(i,0,0)*normfac << std::endl;
        }*/
        ofs.close();
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
    }

    return(EXIT_SUCCESS);
}
