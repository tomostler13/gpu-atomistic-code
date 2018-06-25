// File: bin2Sq2D.cpp
// Author:Tom Ostler
// Created: 25 Feb 2016
// Last-modified: 25 Jun 2018 13:37:43

//The purpose of this section of code is to take a binary file of spin
//files and perform a sum of spins along a given direction (e.g. along z)
//and find the average. Then the 2D fourier transform is performed.

//------------------WORD OF WARNING-----------------------------------
// At the moment this code only works for a cubic lattice and will
// need to be adapted in the future.
//--------------------------------------------------------------------
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
    Array<double> sx,sy,sz,posx,posy,posz;
    Array3D<fftw_complex> s3dx,s3dy,s3dz;
    Array2D<int> ipos;
    Array2D<fftw_complex> Sq2Dx,Sq2Dy,Sq2Dz,Sr2Dx,Sr2Dy,Sr2Dz;
    fftw_plan SrSqx,SrSqy,SrSqz;
    std::cout << "#Reading config file..." << std::flush;
    inputs::readcff(argc,argv);

    double normfac=1./(static_cast<double>(inputs::dimin[0]*inputs::dimin[1]*inputs::dimin[2]));
    std::cout << "Done" << std::endl;

    std::cout << "#Resizing magnetization and spin arrays..." << std::flush;
//    mag.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2],3);
//    mag.IFill(0);
    sx.resize(inputs::nspins);
    sy.resize(inputs::nspins);
    sz.resize(inputs::nspins);
    sx.IFill(0);
    sy.IFill(0);
    sz.IFill(0);
    posx.resize(inputs::nspins);
    posy.resize(inputs::nspins);
    posz.resize(inputs::nspins);
    posx.IFill(0);
    posy.IFill(0);
    posz.IFill(0);
    Sq2Dx.resize(inputs::dimin[0],inputs::dimin[1]);
    Sr2Dx.resize(inputs::dimin[0],inputs::dimin[1]);
    Sr2Dx.IFill(0);Sq2Dx.IFill(0);
    Sq2Dy.resize(inputs::dimin[0],inputs::dimin[1]);
    Sr2Dy.resize(inputs::dimin[0],inputs::dimin[1]);
    Sr2Dy.IFill(0);Sq2Dy.IFill(0);
    Sq2Dz.resize(inputs::dimin[0],inputs::dimin[1]);
    Sr2Dz.resize(inputs::dimin[0],inputs::dimin[1]);
    Sr2Dz.IFill(0);Sq2Dz.IFill(0);
    s3dx.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);s3dx.IFill(0);
    s3dy.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);s3dy.IFill(0);
    s3dz.resize(inputs::dimin[0],inputs::dimin[1],inputs::dimin[2]);s3dz.IFill(0);

    SrSqx=fftw_plan_dft_2d(inputs::dimin[0],inputs::dimin[1],Sr2Dx.ptr(),Sq2Dx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
    SrSqy=fftw_plan_dft_2d(inputs::dimin[0],inputs::dimin[1],Sr2Dy.ptr(),Sq2Dy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
    SrSqz=fftw_plan_dft_2d(inputs::dimin[0],inputs::dimin[1],Sr2Dz.ptr(),Sq2Dz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
    std::cout << "Done" << std::endl;


    //read the positions
    std::ifstream ifspos("pos.bin",std::ios::in | std::ios::binary);
    if(!ifspos.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("could not open position binary file.");
    }
    else
    {
        ifspos.read( (char*)posx.ptr(),inputs::nspins*sizeof(double));
        ifspos.read( (char*)posy.ptr(),inputs::nspins*sizeof(double));
        ifspos.read( (char*)posz.ptr(),inputs::nspins*sizeof(double));
    }
    ifspos.close();
    ipos.resize(inputs::nspins,3);
    ipos.IFill(0);

    //convert to integer grid position
    for(unsigned int i = 0 ; i < inputs::nspins ; i++)
    {
        ipos(i,0)=static_cast<int>(posx[i]/inputs::abc[0]+0.5);
        ipos(i,1)=static_cast<int>(posy[i]/inputs::abc[1]+0.5);
        ipos(i,2)=static_cast<int>(posz[i]/inputs::abc[2]+0.5);
        //std::cout << ipos(i,0) << "\t" << ipos(i,1) << "\t" << ipos(i,2) << std::endl;
    }
    posx.clear();
    posy.clear();
    posz.clear();

    std::cout << "#inputs::dimin = " << inputs::dimin[0] << " , " << inputs::dimin[1] << " , " << inputs::dimin[2] << std::endl;
    for(unsigned int t = inputs::lt ; t < inputs::ut+1 ; t+=inputs::dt)
    {
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        std::stringstream sstr;
        sstr << inputs::fb << t << ".bin";
        std::string str=sstr.str();

        std::cout << "#Opening file " << str << "..." << std::flush;
        std::ifstream ifs(str.c_str(),std::ios::in | std::ios::binary);
        if(!ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open input file");
        }
        std::cout << "done" << std::endl;

        ifs.read( (char*)sx.ptr(),inputs::nspins*sizeof(double));
        ifs.read( (char*)sy.ptr(),inputs::nspins*sizeof(double));
        ifs.read( (char*)sz.ptr(),inputs::nspins*sizeof(double));

        ifs.close();
        if(ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close the binary spin file.");
        }

        std::cout << "#Copying to 3d array..." << std::flush;
        for(int i = 0 ; i < inputs::nspins ; i++)
        {
            int lc[3]={ipos(i,0),ipos(i,1),ipos(i,2)};
            s3dx(lc[0],lc[1],lc[2])[0]=sx[i];
            s3dy(lc[0],lc[1],lc[2])[0]=sy[i];
            s3dz(lc[0],lc[1],lc[2])[0]=sz[i];
        }
        std::cout << "done" << std::endl;


        std::cout << "#Averaging over z..." << std::flush;
        for(int i = 0 ; i < inputs::dimin[0] ; i++)
        {
            for(int j = 0 ; j < inputs::dimin[1] ; j++)
            {
                for(int k = 0 ; k < inputs::dimin[2] ; k++)
                {
                    Sr2Dx(i,j)[0]+=s3dx(i,j,k)[0];
                    Sr2Dy(i,j)[0]+=s3dy(i,j,k)[0];
                    Sr2Dz(i,j)[0]+=s3dz(i,j,k)[0];
                }
                Sr2Dx(i,j)[0]/=static_cast<double>(inputs::dimin[2]);
                Sr2Dy(i,j)[0]/=static_cast<double>(inputs::dimin[2]);
                Sr2Dz(i,j)[0]/=static_cast<double>(inputs::dimin[2]);
            }
        }
        std::cout << "done" << std::endl;
        std::cout << "#Performing FFTs..." << std::flush;
        fftw_execute(SrSqx);
        fftw_execute(SrSqy);
        fftw_execute(SrSqz);
        std::cout << "done" << std::endl;


        std::stringstream osstr;
        osstr << "Sq_" << t << ".dat";
        std::string ostr=osstr.str();
        std::ofstream ofs;
        ofs.open(ostr.c_str());
        if(!ofs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Couldn't open an output file");
        }
        ofs << "#k_x - k_y - Re{Sq_x} - Im{Sq_x} - Re{Sq_y} - Im{Sq_y} - Re{Sq_z} - Im{Sq_z}" << std::endl;

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
                ofs << ri << "\t" << rj << "\t";
                ofs << Sq2Dx(li,lj)[0]*normfac << "\t" << Sq2Dx(li,lj)[1]*normfac << "\t";
                ofs << Sq2Dy(li,lj)[0]*normfac << "\t" << Sq2Dy(li,lj)[1]*normfac << "\t";
                ofs << Sq2Dz(li,lj)[0]*normfac << "\t" << Sq2Dz(li,lj)[1]*normfac << "\t";
                ofs << std::endl;
            }
        }

//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        /*fftw_execute(SxP);
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
        */
        std::cout << "\n" << std::endl;
    }
    std::cout << "-------------------------------------------------------------------------------------------" << std::endl;

    return(EXIT_SUCCESS);
}
