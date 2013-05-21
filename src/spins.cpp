// File: spins.cpp
// Author:Tom Ostler
// Created: 17 Jan 2013
// Last-modified: 21 May 2013 10:47:32
#include <fftw3.h>
#include <libconfig.h++>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/fields.h"
#include "../inc/mat.h"
#include "../inc/intmat.h"
#include "../inc/defines.h"
#include "../inc/maths.h"
#include "../inc/llg.h"
#include "../inc/spins.h"
#include "../inc/util.h"
#include <levmar.h>

#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif


namespace spins
{
    Array3D<fftw_complex> Skx,Sky,Skz;
    Array3D<double> Srx,Sry,Srz;
    Array3D<double> Sznzp;//,Synzp,Sxnzp;
//    Array3D<fftw_complex> Sqznzp,Sqynzp,Sqxnzp;
//    Array3D<double> Szij,Syij,Sxij;
//    Array3D<fftw_complex> Sqzij,Sqyij,Sqxij;
    Array3D<fftw_complex> SpSm;
    unsigned int nzpcplxdim=0;
    double normsize=0;
    double knorm=0;
    double *xdat=NULL;
    util::RunningStat corrLength;
    Array<double> Sx,Sy,Sz,eSx,eSy,eSz;
//    fftw_plan SxP,SyP,SzP,SzcfPF,SycfPF,SxcfPF,SzcfPB,SycfPB,SxcfPB,SzijF,SyijF,SxijF;
    fftw_plan SpSmF;
    unsigned int update=0;
    std::ifstream sfs;
    void initSpins(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Spin details***" << std::endl;
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
        FIXOUT(config::Info,"Resizing arrays:" << std::flush);
        Skx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Sky.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Skz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::cplxdim);
        Srx.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Sry.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Srz.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        Srx.IFill(0);
        Sry.IFill(0);
        Srz.IFill(0);

        Sx.resize(geom::nspins);
        Sy.resize(geom::nspins);
        Sz.resize(geom::nspins);
        SUCCESS(config::Info);
        std::string sc;
        libconfig::Setting &setting = config::cfg.lookup("spins");
        if(llg::ssf)
        {
            FIXOUT(config::Info,"Setting up static structure factor data structures:" <<std::flush);
            SpSm.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
            SpSm.IFill(0);

            SpSmF = fftw_plan_dft_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],SpSm.ptr(),SpSm.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
            std::stringstream sstr2;
            sstr2 << llg::ssffs << ".dat";
            std::string tempstr2=sstr2.str();
            llg::ssffstr.open(tempstr2.c_str());
            if(llg::ssffstr.is_open()!=true)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open file for outputting static structure factor");
            }
            SUCCESS(config::Info);
        }

        if(llg::rscf)
        {
            Sznzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
			Synzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
			Sxnzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
            Szij.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
			Syij.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
			Sxij.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);


            nzpcplxdim=(geom::dim[2]*geom::Nk[2]/2)+1;
            Sqznzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
			Sqynzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
			Sqxnzp.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
            Sqzij.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
			Sqyij.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);
			Sqxij.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],nzpcplxdim);

			SzcfPF = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sznzp.ptr(),Sqznzp.ptr(),FFTW_ESTIMATE);
			SycfPF = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Synzp.ptr(),Sqynzp.ptr(),FFTW_ESTIMATE);
			SxcfPF = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sxnzp.ptr(),Sqxnzp.ptr(),FFTW_ESTIMATE);
			SzijF = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Szij.ptr(),Sqzij.ptr(),FFTW_ESTIMATE);
			SyijF = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Syij.ptr(),Sqyij.ptr(),FFTW_ESTIMATE);
			SxijF = fftw_plan_dft_r2c_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sxij.ptr(),Sqxij.ptr(),FFTW_ESTIMATE);

            SzcfPB = fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sqznzp.ptr(),Sznzp.ptr(),FFTW_ESTIMATE);
            SycfPB = fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sqynzp.ptr(),Synzp.ptr(),FFTW_ESTIMATE);
            SxcfPB = fftw_plan_dft_c2r_3d(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2],Sqxnzp.ptr(),Sxnzp.ptr(),FFTW_ESTIMATE);

            normsize=geom::dim[0]*geom::Nk[0]*geom::dim[0]*geom::Nk[0];
            normsize*=geom::dim[1]*geom::Nk[1]*geom::dim[1]*geom::Nk[1];
            normsize*=geom::dim[2]*geom::Nk[2]*geom::dim[2]*geom::Nk[2];
            knorm=double(geom::nauc)/(double(geom::Nk[0]*geom::Nk[1]*geom::Nk[2]));
            std::stringstream sstr;
            sstr << llg::rscfstr << "all.dat";
            std::string tempstr=sstr.str();
            llg::rscfs.open(tempstr.c_str());
            if(!llg::rscfs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open file for writing correlation function");
            }
            atexit(sexit);
        }


        setting.lookupValue("update",update);
        FIXOUT(config::Info,"Spin update:" << update << " (timesteps)" << std::endl);
        if(update<1)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Spin data should be updated on cpu, i.e. spins::update>0");
        }
        setting.lookupValue("spinconfig",sc);
        FIXOUT(config::Info,"Initial spin config specifier method:" << sc << std::endl);
        if(sc=="file")
        {
            std::string spinfile;
            try
            {
                setting.lookupValue("sf",spinfile);
            }
            catch(const libconfig::SettingTypeException &stex)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Setting type error");
            }
            FIXOUT(config::Info,"Reading spin data from file:" << spinfile << std::flush);
            sfs.open(spinfile.c_str());
            if(!sfs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open file for reading spin data");
            }

            for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                    {
                        int sn=geom::coords(i,j,k,0);
                        if(sn>=0)
                        {
                            int lc[3]={0,0,0};
                            sfs >> lc[0] >> lc[1] >> lc[2];
                            if(lc[0]!=i || lc[1]!=j || lc[2]!=k)
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("The file you are reading in does not correspond to the system dimensions you have specified");
                            }
                            else
                            {
                                sfs >> spins::Sx[sn] >> spins::Sy[sn] >> spins::Sz[sn];
                            }
                        }
                    }
                }
            }
            sfs.close();
            if(sfs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errWarning("Could not close spin file for reading in spin data");
            }
        }
        else if(sc=="align")
        {
            double sr[3]={0,0,0};
            try
            {
                for(unsigned int l = 0 ; l < 3 ; l++){sr[l]=setting["sc"][l];}
            }
            catch(const libconfig::SettingTypeException &stex)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Setting type error");
            }
            FIXOUT(config::Info,"Spins aligned to:" << "(" << sr[0] << "," << sr[1] << "," << sr[2] << ")" <<std::endl);
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                Sx(i)=sr[0];
                Sy(i)=sr[1];
                Sz(i)=sr[2];//sqrt(1.0-(Sx[i]*Sx[i]+Sy[i]*Sy[i]));//+sr[2];
            }
        }
        else if(sc=="random")
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                maths::sphereRandom(Sx(i),Sy(i),Sz(i));
            }
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Spin configuration not recognised, check you config file.");
        }

        FIXOUT(config::Info,"Planning r2c and c2r transforms:" << std::flush);
        //forward transform of spin arrays
        SxP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Srx.ptr(),Skx.ptr(),FFTW_ESTIMATE);
        SyP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Sry.ptr(),Sky.ptr(),FFTW_ESTIMATE);
        SzP=fftw_plan_dft_r2c_3d(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2],Srz.ptr(),Skz.ptr(),FFTW_ESTIMATE);
        SUCCESS(config::Info);

    }
    void FFTForward()
    {
        Skx.IFill(0);
        Sky.IFill(0);
        Skz.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Srx(lc[0],lc[1],lc[2])=Sx[i];
            Sry(lc[0],lc[1],lc[2])=Sy[i];
            Srz(lc[0],lc[1],lc[2])=Sz[i];
            //            std::cout << "Lookup\t" << lc[0] << "\t" << lc[1] << "\t" << lc[2] << "\t" << Skx(lc[0],lc[1],lc[2])[0] << "\t" << Sky(lc[0],lc[1],lc[2])[0] << "\t" << Skz(lc[0],lc[1],lc[2])[0] << std::endl;
        }
        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
    }
    void eFFTForward()
    {
        Skx.IFill(0);
        Sky.IFill(0);
        Skz.IFill(0);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Srx(lc[0],lc[1],lc[2])=eSx[i];
            Sry(lc[0],lc[1],lc[2])=eSy[i];
            Srz(lc[0],lc[1],lc[2])=eSz[i];
            //            std::cout << "Lookup\t" << lc[0] << "\t" << lc[1] << "\t" << lc[2] << "\t" << Skx(lc[0],lc[1],lc[2])[0] << "\t" << Sky(lc[0],lc[1],lc[2])[0] << "\t" << Skz(lc[0],lc[1],lc[2])[0] << std::endl;
        }
        fftw_execute(SxP);
        fftw_execute(SyP);
        fftw_execute(SzP);
    }
    void calcRealSpaceCorrelationFunction(unsigned int t)
    {
//        Sqznzp.IFill(0);
//        Sqynzp.IFill(0);
//        Sqxnzp.IFill(0);
//        Sqzij.IFill(0);
//        Sqyij.IFill(0);
//        Sqxij.IFill(0);
        Sznzp.IFill(0);
//        Synzp.IFill(0);
//        Sxnzp.IFill(0);
//        Szij.IFill(0);
//        Syij.IFill(0);
//        Sxij.IFill(0);

        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            Sznzp(lc[0],lc[1],lc[2])=Sz[i];
            Synzp(lc[0],lc[1],lc[2])=Sy[i];
            Sxnzp(lc[0],lc[1],lc[2])=Sx[i];
            int mc[3]={lc[0],lc[1],lc[2]};
            //if(mc[0]>((geom::dim[0]*geom::Nk[0]/2)-1)){mc[0]=geom::dim[0]*geom::Nk[0]-(mc[0]-1);}
            //if(mc[1]>((geom::dim[1]*geom::Nk[1]/2)-1)){mc[1]=geom::dim[1]*geom::Nk[1]-(mc[1]-1);}
            //if(mc[2]>((geom::dim[2]*geom::Nk[2]/2)-1)){mc[2]=geom::dim[2]*geom::Nk[2]-(mc[2]-1);}
//            std::cout << lc[0] << "\t" << lc[1] << "\t" << lc[2] << "\t" << mc[0] << "\t" << mc[1] << "\t" << mc[2] << std::endl;std::cin.get();
            Szij(mc[0],mc[1],mc[2])=Sz[i];
            Syij(mc[0],mc[1],mc[2])=Sy[i];
            Sxij(mc[0],mc[1],mc[2])=Sx[i];
        }
        for( int i = 0 ; i < geom::dim[0]*geom::Nk[0]/2 ; i++)
        {
            for( int j = 0 ; j < geom::dim[1]*geom::Nk[1]/2 ; j++)
            {
                for( int k = 0 ; k < geom::dim[2]*geom::Nk[2]/2 ; k++)
                {
                    for( int x = 0 ; x < geom::dim[0]*geom::Nk[0]/2 ; x++)
                    {
                        for( int y = 0 ; y < geom::dim[1]*geom::Nk[1]/2 ; y++)
                        {
                            for( int z = 0 ; z < geom::dim[2]*geom::Nk[2]/2 ; z++)
                            {
                                int sn=geom::coords(x,y,z,0);
                                int rij[3]={i-x,j-y,k-z};
                                int nrij[3]={0,0,0};
                                for(unsigned int c = 0 ; c < 3 ; c++)
                                {
                                    if(rij[c]<0)//geom::dim[c]*geom::Nk[c]/2)
                                    {
                                        nrij[c]=geom::dim[c]*geom::Nk[c]+rij[c];
                                    }
                                    Sxij(nrij[0],nrij[1],nrij[2])=Sx[sn];
                                    Syij(nrij[0],nrij[1],nrij[2])=Sy[sn];
                                    Szij(nrij[0],nrij[1],nrij[2])=Sz[sn];
                                }
                            }
                        }
                    }
                }
            }
        }
        fftw_execute(SzcfPF);
        fftw_execute(SycfPF);
        fftw_execute(SxcfPF);
        fftw_execute(SxijF);
        fftw_execute(SyijF);
        fftw_execute(SzijF);


        for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < nzpcplxdim ; k++)
                {
                    const double Sqz[2]={Sqznzp(i,j,k)[0],Sqznzp(i,j,k)[1]};
                    const double Sqy[2]={Sqynzp(i,j,k)[0],Sqynzp(i,j,k)[1]};
                    const double Sqx[2]={Sqxnzp(i,j,k)[0],Sqxnzp(i,j,k)[1]};
                    const double CCSqz[2]={Sqzij(i,j,k)[0],Sqzij(i,j,k)[1]};
					const double CCSqy[2]={Sqyij(i,j,k)[0],Sqyij(i,j,k)[1]};
					const double CCSqx[2]={Sqxij(i,j,k)[0],Sqxij(i,j,k)[1]};
                    Sqznzp(i,j,k)[0]=(Sqx[0]*CCSqx[0]-Sqx[1]*CCSqx[1]);
                    Sqynzp(i,j,k)[0]=(Sqy[0]*CCSqy[0]-Sqy[1]*CCSqy[1]);
                    Sqxnzp(i,j,k)[0]=(Sqz[0]*CCSqz[0]-Sqz[1]*CCSqz[1]);//Sqznzp(i,j,k)[0]*Sqznzp(i,j,k)[0]-Sqznzp(i,j,k)[1]*Sqznzp(i,j,k)[1];
                    Sqznzp(i,j,k)[1]=(Sqx[0]*CCSqx[1]+Sqx[1]*CCSqx[0]);
                    Sqynzp(i,j,k)[1]=(Sqy[0]*CCSqy[1]+Sqy[1]*CCSqy[0]);
                    Sqxnzp(i,j,k)[1]=(Sqz[0]*CCSqz[1]+Sqz[1]*CCSqz[0]);//Sqznzp(i,j,k)[0]*Sqznzp(i,j,k)[1]+Sqznzp(i,j,k)[1]*Sqznzp
                }
            }
        }
        fftw_execute(SzcfPB);
        fftw_execute(SycfPB);
        fftw_execute(SxcfPB);
        /* CAN BE USED TO OUTPUT THE WHOLE CORRELATION FUNCTION IN ALL X,Y,Z directions
           for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
           {
           for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
           {
           for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
           {
           if(geom::coords(i,j,k,0)>-1)
           {
           int ijk[3]={i,j,k};
           for(unsigned int c = 0 ; c < 3 ; c++)
           {
           if(ijk[c]>geom::dim[c]*geom::Nk[c]/2)
           {
           ijk[c]=ijk[c]-geom::dim[c]*geom::Nk[c];
           }
           }
           llg::rscfs << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << "\t" << Sznzp(i,j,k)/(normsize*knorm) << std::endl;
           }
           }
           }
           llg::rscfs << std::endl;
           }
           std::cin.get();*/
        for(int k = -(geom::dim[2]*geom::Nk[2]/2)+1 ; k <= 0 ; k++)
        {

            int arluv=geom::dim[2]*geom::Nk[2]+k;
            if(geom::coords(0,0,arluv,0)>-1)
            {
                llg::rscfs << t << "\t"<< k << "\t" << (Sznzp(0,0,arluv)+Sxnzp(0,0,arluv)+Synzp(0,0,arluv))/(normsize*knorm) << std::endl;
            }
        }
        for(int k = 1 ; k < (geom::dim[2]*geom::Nk[2]/2) ; k++)
        {
            if(geom::coords(0,0,k,0)>-1)
            {
                llg::rscfs << t << "\t" << k << "\t" << (Sznzp(0,0,k)+Synzp(0,0,k)+Sxnzp(0,0,k))/(normsize*knorm) << std::endl;
            }
        }
        llg::rscfs << std::endl;// << std::endl;
        //fitcorr(t);
        //llg::rscfs.close();
        //if(llg::rscfs.is_open())
        //{
        //    error::errPreamble(__FILE__,__LINE__);
        //    error::errWarning("Could not close file for writing correlation function");
        //}
    }
    void calcStaticStructureFactor(unsigned int t)
    {
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int lc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            SpSm(lc[0],lc[1],lc[2])[0]=Sx[i];
            SpSm(lc[0],lc[1],lc[2])[1]=Sy[i];
        }
        fftw_execute(SpSmF);
        if(geom::dim[2]>1)
        {
            //X -> Gamma within the brillouin zone. We are only currently interested in
            //positive k-vectors
            for(int k = (geom::dim[2]*geom::Nk[2]/2)-1 ; k > 0 ; k--)
            {
                llg::ssffstr << t << "\t" << -k << "\t" << SpSm(0,0,k)[0] << "\t" << SpSm(0,0,k)[1] << std::endl;
            }
        }
        if(geom::dim[0] == geom::dim[2] && geom::dim[1]==geom::dim[2])
        {
            //Gamma -> M
            for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2]/2 ; k++)
            {
                llg::ssffstr << t << "\t" << k << "\t" << SpSm(0,k,k)[0] << "\t" << SpSm(0,k,k)[1] << std::endl;
            }
        }
        llg::ssffstr << std::endl;
    }
    void corrfunc(double *p,double *ydat,int m,int n,void *data)
    {
        register int i;
        for(int i = 0 ; i < n ; i++)
        {
            ydat[i]=udf(p,xdat[i]);//p[0]*exp(-double(i)/p[1])+p[2];
        }
    }
    double udf(double *p,double x)
    {
        return(p[0]*exp(-(x-p[1])/p[2])+p[3]);
    }

    struct xtradata{
        char msg[128];
    };

    void fitcorr(unsigned int t)
    {
        const int m=4,ndp=((geom::dim[2]*geom::Nk[2]/2)-geom::Nk[2])/geom::Nk[2];
        double ydat[ndp],p[m],opts[LM_OPTS_SZ],info[LM_INFO_SZ];
        int ret;
        struct xtradata data;
        unsigned int count=0;
        xdat=new double [ndp];
        int initk=0;
        bool ikc=false;
        for(int k = 1 ; k < (geom::dim[2]*geom::Nk[2]/2) ; k++)
        {
            if(geom::coords(0,0,k,0)>-1)
            {
                if(ikc==false)
                {
                    initk=k;
                    ikc=true;
                }
                xdat[count]=k-initk;
                ydat[count]=Sznzp(0,0,k)/(normsize*knorm);
                count++;
            }
        }
        if(count!=ndp)
        {
            std::cout << "Number of data points = " << count << "\tshould be " << ndp << std::endl;
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Number of data points not correct");
        }
        //initial guess for parameters
        p[0]=Sznzp(0,0,initk);
        p[1]=0.00001;//2.0*M_PI/(4.0*geom::dim[2]*geom::Nk[2]);
        p[2]=geom::dim[2]*geom::Nk[2]/2.0;
        p[3]=0.001;
        /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
        opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
        opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used
        ret=dlevmar_dif(corrfunc,p,ydat,m,ndp,1000,opts,info,NULL,NULL,(void *)&data);
        //std::cout << t << "\t" << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << p[3] << std::endl;
        if(p[2]<geom::dim[2]*geom::Nk[2] && p[2] > 0)
        {
            corrLength.Push(p[2]);
        }
        //std::cout << ret << std::endl;
        //std::cerr << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << p[3] << std::endl;
        ///* CAN BE USED FOR DEBUGGING FITTING.
        /*double lv=0.0;//double(initk);
          double dv=(geom::dim[2]*geom::Nk[2]/2-geom::Nk[2])/1000.0;
          for(unsigned int i = 0 ; i < 1000 ; i++)
          {
          std::cerr << lv+double(i)*dv << "\t" << udf(p,lv+double(i)*dv) << std::endl;
          }
          std::cin.get();*/

        delete [] xdat;
        xdat=NULL;

        //printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
        //printf("Best fit parameters: %.7g %.7g %.7g %.7g\n", p[0], p[1], p[2], p[3]);

    }

    void sexit(void)
    {
        llg::rscfs.close();
        if(llg::rscfs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close file for writing correlation function");
        }
        FIXOUT(config::Info,"Average correlation length (mean , variance):" << spins::corrLength.Mean() << "\t" << spins::corrLength.Variance() << std::endl);
        std::cout << "#Average correlation length = " << spins::corrLength.Mean() << std::endl;
    }


}
