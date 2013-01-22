// File: util.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 22 Jan 2013 11:53:09
// Contains useful functions and classes
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/intmat.h"
#include "../inc/geom.h"
#include <omp.h>
extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}
namespace util
{
    void inverse(double* A, int N)
    {
        int *IPIV = new int[N+1];
        int LWORK = N*N;
        double *WORK = new double[LWORK];
        int INFO;

        dgetrf_(&N,&N,A,&N,IPIV,&INFO);
        dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

        delete IPIV;
        delete WORK;
    }
    void cpuConvFourier()
    {
        fields::Hkx.IFill(0);
        fields::Hky.IFill(0);
        fields::Hkz.IFill(0);

        //perform convolution in fourier space
        for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::cplxdim ; k++)
                {

                    fields::Hkx(i,j,k)[0]=intmat::Nxx(i,j,k)[0]*spins::Skx(i,j,k)[0]-intmat::Nxx(i,j,k)[1]*spins::Skx(i,j,k)[1]
                        +intmat::Nxy(i,j,k)[0]*spins::Sky(i,j,k)[0]-intmat::Nxy(i,j,k)[1]*spins::Sky(i,j,k)[1]
                        +intmat::Nxz(i,j,k)[0]*spins::Skz(i,j,k)[0]-intmat::Nxz(i,j,k)[1]*spins::Skz(i,j,k)[1];
                    fields::Hkx(i,j,k)[1]=intmat::Nxx(i,j,k)[0]*spins::Skx(i,j,k)[1]+intmat::Nxx(i,j,k)[1]*spins::Skx(i,j,k)[0]
                        +intmat::Nxy(i,j,k)[0]*spins::Sky(i,j,k)[1]+intmat::Nxy(i,j,k)[1]*spins::Sky(i,j,k)[0]
                        +intmat::Nxz(i,j,k)[0]*spins::Skz(i,j,k)[1]+intmat::Nxz(i,j,k)[1]*spins::Skz(i,j,k)[0];

                    fields::Hky(i,j,k)[0]=intmat::Nyx(i,j,k)[0]*spins::Skx(i,j,k)[0]-intmat::Nyx(i,j,k)[1]*spins::Skx(i,j,k)[1]
                        +intmat::Nyy(i,j,k)[0]*spins::Sky(i,j,k)[0]-intmat::Nyy(i,j,k)[1]*spins::Sky(i,j,k)[1]
                        +intmat::Nyz(i,j,k)[0]*spins::Skz(i,j,k)[0]-intmat::Nyz(i,j,k)[1]*spins::Skz(i,j,k)[1];
                    fields::Hky(i,j,k)[1]=intmat::Nyx(i,j,k)[0]*spins::Skx(i,j,k)[1]+intmat::Nyx(i,j,k)[1]*spins::Skx(i,j,k)[0]
                        +intmat::Nyy(i,j,k)[0]*spins::Sky(i,j,k)[1]+intmat::Nyy(i,j,k)[1]*spins::Sky(i,j,k)[0]
                        +intmat::Nyz(i,j,k)[0]*spins::Skz(i,j,k)[1]+intmat::Nyz(i,j,k)[1]*spins::Skz(i,j,k)[0];

                    fields::Hkz(i,j,k)[0]=intmat::Nzx(i,j,k)[0]*spins::Skx(i,j,k)[0]-intmat::Nzx(i,j,k)[1]*spins::Skx(i,j,k)[1]
                        +intmat::Nzy(i,j,k)[0]*spins::Sky(i,j,k)[0]-intmat::Nzy(i,j,k)[1]*spins::Sky(i,j,k)[1]
                        +intmat::Nzz(i,j,k)[0]*spins::Skz(i,j,k)[0]-intmat::Nzz(i,j,k)[1]*spins::Skz(i,j,k)[1];
                    fields::Hkz(i,j,k)[1]=intmat::Nzx(i,j,k)[0]*spins::Skx(i,j,k)[1]+intmat::Nzx(i,j,k)[1]*spins::Skx(i,j,k)[0]
                        +intmat::Nzy(i,j,k)[0]*spins::Sky(i,j,k)[1]+intmat::Nzy(i,j,k)[1]*spins::Sky(i,j,k)[0]
                        +intmat::Nzz(i,j,k)[0]*spins::Skz(i,j,k)[1]+intmat::Nzz(i,j,k)[1]*spins::Skz(i,j,k)[0];

                }
            }
        }
    }
}

