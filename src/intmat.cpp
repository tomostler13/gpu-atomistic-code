// File: intmat.cpp
// Author:Tom Ostler
// Created: 16 Jan 2012
// Last-modified: 18 Oct 2016 12:52:56
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <libconfig.h++>
#include "../inc/intmat.h"
#include "../inc/arrays.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/unitcell.h"
#include "../inc/llg.h"
#include "../inc/defines.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace intmat
{
    Array7D<fftw_complex> Nkab;
    Array7D<fftw_complex> Nrab;
    Array5D<fftw_complex> dipNkab,dipNrab;
    //These are for the hybrid method. Each plane in these matrices
    //stores the set of interaction matrices for the planes
    Array5D<fftw_complex> hNkab,hNrab;
    Array<unsigned int> zpsn;
    fftw_plan ftP,dipftP;

    void initIntmat()
    {
        config::Info.width(45);config::Info << std::right << "*" << "**Interaction matrix details***" << std::endl;
        if(geom::Nkset==false)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Cannot have fft method without specifying Nm");
        }
        //For the interaction matrices we can use a 7D array to make the loops very compact
        //The first element is the species (i)
        //Second is the species (j) that species (i) is interacting with
        //Third and fourth elements are the elements of the tensor
        //fifth, sixth and seventh are the space (reciprocal space) elements
        if(config::exchm<99)//the use the FFT for everything
        {
            Nkab.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),3,3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            Nrab.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),3,3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            FIXOUT(config::Info,"The real space interaction matrix contains:" << Nrab.size() << " elements (double complex)" << std::endl);
            FIXOUT(config::Info,"The reciprocal space interaction matrix contains:" << Nkab.size() << " elements (double complex)" << std::endl);
            Nrab.IFill(0.0);
            Nkab.IFill(0.0);
            config::Info.width(45);config::Info << std::right << "*" << "**Interaction matrix details***" << std::endl;
            FIXOUT(config::Info,"Setting up fftw of interaction matrix:" << std::endl);
            //FOR DEBUGGING: TO LOOK AT THE REAL SPACE INTERACTION MATRIX
            /*for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
              {
              for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
              {
              for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
              {
              std::cout << i << "\t" << j << "\t" << k << "\t" << Nrab(0,0,0,0,i,j,k)[0] << "\t"<< Nrab(0,0,1,1,i,j,k)[0] << "\t"<< Nrab(0,0,2,2,i,j,k)[0] << std::endl;
              }
              }
              }
              std::cin.get();*/
            //set the time limit for optimizing the plan of the FFT to 300 seconds
            fftw_set_timelimit(300);
            int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zps;
            int idist=geom::zps;

            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of interaction matrix" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3 << std::endl);
            FIXOUT(config::Log,"Pointer of real space int mat (Nrab):" << Nrab.ptr() << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUT(config::Log,"Pointer of reciprocal space int mat (Nkab):" << Nkab.ptr() << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) = " << "FFTW_FORWARD" << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            int howmany=geom::ucm.GetNMS()*geom::ucm.GetNMS()*3*3;
            ftP=fftw_plan_many_dft(3,n,howmany,Nrab.ptr(),inembed,istride,idist,Nkab.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_PATIENT);
        }
        else if(config::exchm>98)//here we are using the hybrid scheme
        {
            hNrab.resize(3,3,geom::dim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            hNkab.resize(3,3,geom::dim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
            hNrab.IFill(0);
            hNkab.IFill(0);
            config::Info.width(45);config::Info << std::right << "*" << "**Interaction matrix details***" << std::endl;
            FIXOUT(config::Info,"Setting up fftw of set of 2D interaction matrices (see log):" << std::flush);
            fftw_set_timelimit(300);
            int n[2]={geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
            int *inembed=n;
            int *onembed=n;
            int istride=1;
            int ostride=1;
            int odist=geom::zpdim[1]*geom::Nk[1]*geom::zpdim[2]*geom::Nk[2];
            int idist=odist;
            config::openLogFile();
            config::printline(config::Log);
            FIXOUT(config::Log,"Parameters entering into FFTW plan of interaction matrix" << std::endl);
            FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],0);
            FIXOUT(config::Log,"rank (dimension of FFT) = " << 2 << std::endl);
            FIXOUT(config::Log,"How many (FFT's) = " << geom::Nk[0]*geom::dim[0]*3*3 << std::endl);
            FIXOUT(config::Log,"Pointer of real space int mat (hNrab):" << hNrab.ptr() << std::endl);
            FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],0);
            FIXOUT(config::Log,"istride = " << istride << std::endl);
            FIXOUT(config::Log,"idist = " << idist << std::endl);
            FIXOUT(config::Log,"Pointer of reciprocal space int mat (hNkab):" << hNkab.ptr() << std::endl);
            FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],0);
            FIXOUT(config::Log,"ostride = " << ostride << std::endl);
            FIXOUT(config::Log,"odist = " << odist << std::endl);
            FIXOUT(config::Log,"Direction (sign) = " << "FFTW_FORWARD" << std::endl);
            FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
            int howmany=geom::Nk[0]*geom::dim[0]*3*3;
            ftP=fftw_plan_many_dft(2,n,howmany,hNrab.ptr(),inembed,istride,idist,hNkab.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_PATIENT);
            SUCCESS(config::Info);

        }


    }
    void initDipIntmat()
    {
        config::Info.width(45);config::Info << std::right << "*" << "**Interaction matrix details***" << std::endl;
        dipNkab.resize(3,3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        dipNrab.resize(3,3,geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
        FIXOUT(config::Info,"The real space interaction matrix contains:" << dipNrab.size() << " elements (double complex)" << std::endl);
        FIXOUT(config::Info,"The reciprocal space interaction matrix contains:" << dipNkab.size() << " elements (double complex)" << std::endl);
        dipNrab.IFill(0.0);
        dipNkab.IFill(0.0);
        config::Info.width(45);config::Info << std::right << "*" << "**Interaction matrix details***" << std::endl;
        FIXOUT(config::Info,"Setting up fftw of interaction matrix:" << std::endl);
        //FOR DEBUGGING: TO LOOK AT THE REAL SPACE INTERACTION MATRIX
        /*for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    std::cout << i << "\t" << j << "\t" << k << "\t" << Nrab(0,0,0,0,i,j,k)[0] << "\t"<< Nrab(0,0,1,1,i,j,k)[0] << "\t"<< Nrab(0,0,2,2,i,j,k)[0] << std::endl;
                }
            }
        }
        std::cin.get();*/
        //set the time limit for optimizing the plan of the FFT to 60 seconds
        fftw_set_timelimit(60);
        int n[3]={geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]};
        int *inembed=n;
        int *onembed=n;
        int istride=1;
        int ostride=1;
        int odist=geom::zps;
        int idist=geom::zps;

        config::openLogFile();
        config::printline(config::Log);
        FIXOUT(config::Log,"Parameters entering into FFTW plan of interaction matrix" << std::endl);
        FIXOUTVEC(config::Log,"Dimensions of FFT = ",n[0],n[1],n[2]);
        FIXOUT(config::Log,"rank (dimension of FFT) = " << 3 << std::endl);
        FIXOUT(config::Log,"How many (FFT's) = " << 3*3 << std::endl);
        FIXOUT(config::Log,"Pointer of real space int mat (dipNrab):" << dipNrab.ptr() << std::endl);
        FIXOUTVEC(config::Log,"inembed = ",inembed[0],inembed[1],inembed[2]);
        FIXOUT(config::Log,"istride = " << istride << std::endl);
        FIXOUT(config::Log,"idist = " << idist << std::endl);
        FIXOUT(config::Log,"Pointer of reciprocal space int mat (dipNkab):" << dipNkab.ptr() << std::endl);
        FIXOUTVEC(config::Log,"onembed = ",onembed[0],onembed[1],onembed[2]);
        FIXOUT(config::Log,"ostride = " << ostride << std::endl);
        FIXOUT(config::Log,"odist = " << odist << std::endl);
        FIXOUT(config::Log,"Direction (sign) = " << "FFTW_FORWARD" << std::endl);
        FIXOUT(config::Log,"flags = " << "FFTW_PATIENT" << std::endl);
        int howmany=3*3;
        dipftP=fftw_plan_many_dft(3,n,howmany,dipNrab.ptr(),inembed,istride,idist,dipNkab.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_PATIENT);


    }

    void fillIntmat()
    {
        FIXOUT(config::Info,"Filling interaction matrix with dipole-dipole terms:" << std::flush);
        //atom counter
        unsigned ac=0;
        //loop over the unit cells
        double rij[3]={0,0,0};
        //the identity matrix
        int I[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        //check if we have a single layer of atoms in any dimension
        //we are not going to output any information about this at this
        //point as we are going to do it in the exchange
        bool checkmonolayer[3]={false,false,false};
        for(unsigned int xyz = 0 ; xyz < 3; xyz++)
        {
            if(geom::dim[xyz]*geom::Nk[xyz] < 2)
            {
                checkmonolayer[xyz]=true;
            }
        }
        //loop over species 1
        for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
        {
            for(unsigned int s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
            {
                for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
                {
                    for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
                    {
                        for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                        {
                            if(!(i==0 && j==0 && k==0))
                            {
                                int lc[3]={i,j,k};
                                int tc[3]={lc[0],lc[1],lc[2]};
                                //check if we have a single layer on any dimension
                                if((abs(tc[0])>0 && checkmonolayer[0]==true) || (abs(tc[1])>0 && checkmonolayer[1]==true) || (abs(tc[2])>0 && checkmonolayer[2]==true) )
                                {
                                    //don't add anything to the interaction matrix
                                }
                                else
                                {
                                    //The interaction matrix must be wrapped around for a C-array format
                                    for(unsigned int l = 0 ; l < 3 ; l++)
                                    {
                                        if(static_cast<unsigned int>(lc[l])>geom::dim[l]*geom::Nk[l])
                                        {
                                            lc[l]=geom::dim[l]*geom::Nk[l]-lc[l];
                                            tc[l]=geom::zpdim[l]*geom::Nk[l]+lc[l];
                                        }
                                        rij[l]=static_cast<double>(lc[l])*geom::abc[l]/static_cast<double>(geom::Nk[l]);
                                    }
                                    double mrij=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                                    double oomrij3=1./(mrij*mrij*mrij);
                                    //unit vector from i to j
                                    double eij[3]={rij[0]/mrij,rij[1]/mrij,rij[2]/mrij};
                                    //loop over (alpha) the row of the tensor
                                    for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                                    {
                                        //loop over the column of the tensor (alpha)
                                        for(unsigned int beta = 0 ; beta < 3 ; beta++)
                                        {
                                            Nrab(s1,s2,alpha,beta,tc[0],tc[1],tc[2])[0]+=1e-7*((3.0*eij[alpha]*eij[beta])-I[alpha][beta])*oomrij3*geom::ucm.GetMuBase(s2)*llg::muB;
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }
        config::Info << "Done" << std::endl;
    }
    void fillDipIntmat()
    {
        FIXOUT(config::Info,"Filling interaction matrix with dipole-dipole terms:" << std::flush);
        //atom counter
        unsigned ac=0;
        //loop over the unit cells
        double rij[3]={0,0,0};
        //the identity matrix
        int I[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        //check if we have a single layer of atoms in any dimension
        //we are not going to output any information about this at this
        //point as we are going to do it in the exchange
        bool checkmonolayer[3]={false,false,false};
        for(unsigned int xyz = 0 ; xyz < 3; xyz++)
        {
            if(geom::dim[xyz]*geom::Nk[xyz] < 2)
            {
                checkmonolayer[xyz]=true;
            }
        }
        for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    if(!(i==0 && j==0 && k==0))
                    {
                        int lc[3]={i,j,k};
                        int tc[3]={lc[0],lc[1],lc[2]};
                        //check if we have a single layer on any dimension
                        if((abs(tc[0])>0 && checkmonolayer[0]==true) || (abs(tc[1])>0 && checkmonolayer[1]==true) || (abs(tc[2])>0 && checkmonolayer[2]==true) )
                        {
                            //don't add anything to the interaction matrix
                        }
                        else
                        {
                            //The interaction matrix must be wrapped around for a C-array format
                            for(unsigned int l = 0 ; l < 3 ; l++)
                            {
                                if(static_cast<unsigned int>(lc[l])>geom::dim[l]*geom::Nk[l])
                                {
                                    lc[l]=geom::dim[l]*geom::Nk[l]-lc[l];
                                    tc[l]=geom::zpdim[l]*geom::Nk[l]+lc[l];
                                }
                                rij[l]=static_cast<double>(lc[l])*geom::abc[l]/static_cast<double>(geom::Nk[l]);
                            }
                            double mrij=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                            double oomrij3=1./(mrij*mrij*mrij);
                            //unit vector from i to j
                            double eij[3]={rij[0]/mrij,rij[1]/mrij,rij[2]/mrij};
                            //loop over (alpha) the row of the tensor
                            for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                            {
                                //loop over the column of the tensor (alpha)
                                for(unsigned int beta = 0 ; beta < 3 ; beta++)
                                {
                                    //Here we do not multiply by the moment of the neigbour, we incorporate that into the spin array
                                    dipNrab(alpha,beta,tc[0],tc[1],tc[2])[0]+=1e-7*((3.0*eij[alpha]*eij[beta])-I[alpha][beta])*oomrij3*llg::muB;
                                }
                            }
                        }
                    }
                }

            }
        }
        config::Info << "Done" << std::endl;
    }
    void fftIntmat()
    {
        config::printline(config::Info);
        FIXOUT(config::Info,"Performing forward transform of N and destroying plan:" << std::flush);
        //execute the demag tensor transforms and destroy the plans
        fftw_execute(ftP);
        /*for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    std::cout << i << "\t" << j << "\t" << k << "\t" << Nkab(0,0,0,0,i,j,k)[0] << "\t"<< Nkab(0,0,1,1,i,j,k)[0] << "\t"<< Nkab(0,0,2,2,i,j,k)[0] << std::endl;
                }
            }
        }
        std::cin.get();*/
        if(config::exchm<99)
        {
            Nrab.clear();
        }
        else
        {
            for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                    {
                        //std::cout << i << "\t" << j << "\t" << k << "\t" << hNkab(0,0,i,j,k)[0] << "\t"<< hNkab(1,1,i,j,k)[0] << "\t"<< hNkab(2,2,i,j,k)[0] << std::endl;
                    }
                }
            }
            //std::cin.get();
            hNrab.clear();
        }
        //destroy the plans
        fftw_destroy_plan(ftP);
        config::Info << "Done" << std::endl;

    }
    void fftDipIntmat()
    {
        config::printline(config::Info);
        FIXOUT(config::Info,"Performing forward transform of N and destroying plan:" << std::flush);
        //execute the demag tensor transforms and destroy the plans
        fftw_execute(dipftP);
        /*for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    std::cout << i << "\t" << j << "\t" << k << "\t" << Nkab(0,0,0,0,i,j,k)[0] << "\t"<< Nkab(0,0,1,1,i,j,k)[0] << "\t"<< Nkab(0,0,2,2,i,j,k)[0] << std::endl;
                }
            }
        }
        std::cin.get();*/
        dipNrab.clear();
        //destroy the plans
        fftw_destroy_plan(dipftP);
        config::Info << "Done" << std::endl;

    }


}
