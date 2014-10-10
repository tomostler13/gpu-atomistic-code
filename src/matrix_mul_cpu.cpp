//File matrix_mul_cpu.cpp
// Author: Tom Ostler
// Created: 08 Oct 2014
// Last-modified: 10 Oct 2014 14:57:12
#include "../inc/arrays.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include <string>
#include <iostream>
#include <cstdlib>
namespace matmul
{
    void spmv_dia_diag(unsigned int uN,int num_diags,Array<int>& offset,Array<double>& dataxx,Array<double>& datayy,Array<double>& datazz,Array<double>& Sx,Array<double>& Sy,Array<double>& Sz,Array<double>& Hx,Array<double>& Hy,Array<double>& Hz)
    {
        int N=static_cast<int>(uN);
        //loop over nspins
        for(int i = 0 ; i < static_cast<int>(N) ; i++)
        {
//            Hx[i]=0;
//            Hy[i]=0;
//            Hz[i]=0;
std::cout << "spin i = " << i <<std::endl;
            double sumx=0.0,sumy=0.0,sumz=0.0;
            for(int n = 0 ; n < num_diags ; n++)
            {
                int j=i+offset[n];
                //int arlu=i*num_diags+n;
                int arlu=N*n+i;

                double valx=dataxx[arlu],valy=datayy[arlu],valz=datazz[arlu];
                if(j >= 0 && j < N)
                {

                    std::cout << "neighbour " << j << "\t" << datazz[arlu] << std::endl;
                    sumx+=(valx*Sx[j]);
                    sumy+=(valy*Sy[j]);
                    sumz+=(valz*Sz[j]);
                }
            }
            //multiply by 1/2 to account for double counting
            Hx[i]+=sumx;
            Hy[i]+=sumy;
            Hz[i]+=sumz;

        }
    }
    void spmv_dia_offdiag(unsigned int uN,int num_diags,Array<int>& offset,Array<double>& dataxy,Array<double>& dataxz,Array<double>& datayx,Array<double>& datayz,Array<double>& datazx,Array<double>& datazy,Array<double>& Sx,Array<double>& Sy,Array<double>& Sz,Array<double>& Hx,Array<double>& Hy,Array<double>& Hz)
    {
        int N=static_cast<int>(uN);
        //loop over nspins
/*        for(int i = 0 ; i < static_cast<int>(N) ; i++)
        {
            double sumx=0.0,sumy=0.0,sumz=0.0;
            for(int n = 0 ; n < num_diags ; n++)
            {
                int j=i+offset[n];
                //int arlu=i*num_diags+n;
                int arlu=N*n+i;
                double valx=dataxx[arlu],valy=datayy[arlu],valz=datazz[arlu];
                if(j >= 0 && j < N)
                {

                    sumx+=(valx*Sx[j]);
                    sumy+=(valy*Sy[j]);
                    sumz+=(valz*Sz[j]);
                }
            }
            //multiply by 1/2 to account for double counting
            Hx[i]+=sumx;
            Hy[i]+=sumy;
            Hz[i]+=sumz;

        }
        */
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The anti-symmetric exchange is not currently coded for use with the DIA sparse matrix multiplication method.");
    }

    void spmv_csr_diag(unsigned int N,Array<unsigned int>& xadj,Array<unsigned int>& adjncy,Array<double>& dataxx,Array<double>& datayy,Array<double>& datazz,Array<double>& Sx,Array<double>& Sy,Array<double>& Sz,Array<double>& Hx,Array<double>& Hy,Array<double>& Hz )
    {
        for(unsigned int i = 0 ; i < N ; i++)
        {
            double sumx=0.0,sumy=0.0,sumz=0.0;
            for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
            {
                unsigned int neigh=adjncy[j];
                sumx+=dataxx[j]*Sx[j];
                sumy+=datayy[j]*Sy[j];
                sumz+=datazz[j]*Sz[j];
            }
            Hx[i]+=sumx;
            Hy[i]+=sumy;
            Hz[i]+=sumz;
        }
    }

    void spmv_csr_offdiag()
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The anti-symmetric exchange is not currently coded for use with the CSR sparse matrix multiplication method.");
    }

}
