//File matrix_conv.cpp
// Author: Tom Ostler
// Created: 07 Jan 2014
// Last-modified: 07 Oct 2014 13:40:49
// The routines within this file convert a 2D matrix to
// a number of formats depending on the routine used. The
// return structures depend on the storage format
#include "../inc/arrays.h"
#include "../inc/matrix_conv.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include <string>
#include <iostream>
#include <cstdlib>
namespace matconv
{

    // This calculate the diagonal offsets for the DIA format
    // and determine the data.
    // The first arguement is the array to store the offsets
    // The second is the J matrix (4 dimensions because it
    // represents all elements of the interactio matrix).
    // This function is overloaded. If the arguement list
    // contains only one offset array then we calculate that
    // for the diagonals of the exchange tensor.
    void conv_intmat_to_dia(Array<int>& offset,Array4D<double>& JMat,unsigned int& num_diag,unsigned int& nN,Array<double>& data)
    {
        //for counting the number of non-zero diagonals
        unsigned int dia_count=0;
        int N=static_cast<int>(nN);
        //set the number of diagonals to zero
        std::cout << "N=" << N << std::endl;
        num_diag=0;
        for(int i = -N+1 ; i < N ; i++)
        {
            //calculate the number of elements in the diagonal
            int num_elem=0;
            if(i<=0)
            {
                num_elem=i+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-i;
            }
            unsigned int istart=0,jstart=0;
            if(i<=0)
            {
                istart=-i;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=i;
            }
            int iloop=istart,jloop=jstart;

            double abs_J=0.0;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {
                abs_J+=fabs(JMat(iloop,jloop,0,0));
                iloop++;
                jloop++;
            }
            if(abs_J>1e-30)
            {
                num_diag++;
            }
        }
        //resize the offset array
        offset.resize(num_diag);
        data.resize(num_diag*N);
        data.IFill(0);
        offset.IFill(0);
        unsigned int nd=0;
        //loop over the diagonals again and set the value of the offset
        for(int i = -N+1 ; i < N ; i++)
        {
            //calculate the number of elements in the diagonal
            int num_elem=0,num_zeros=0;
            if(i<=0)
            {
                num_elem=i+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-i;
            }
            num_zeros=N-num_elem;
            unsigned int istart=0,jstart=0;
            if(i<=0)
            {
                istart=-i;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=i;
            }
            int iloop=istart,jloop=jstart;


            double abs_J=0.0;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {

                abs_J+=fabs(JMat(iloop,jloop,0,0));
                iloop++;
                jloop++;
            }
            if(abs_J>1e-30)
            {
                iloop=istart;
                jloop=jstart;
                for(unsigned int j = 0 ; j < num_elem ; j++)
                {
                    if(i<0)
                    {
                        data[num_diag*nd+(num_zeros)+j]=JMat(iloop,jloop,0,0);
                    }
                    else
                    {
                        data[num_diag*nd+j]=JMat(iloop,jloop,0,0);
                    }
                    iloop++;
                    jloop++;
                }
                offset[nd]=i;
                nd++;
            }
        }
        //for debugging the data array
        /*
        for(unsigned int i = 0 ; i < num_diag ; i++)
        {
            int os=offset[i];
            std::cout << "Offset " << i << " is " << os << std::endl;
            for(unsigned int j = 0 ; j < N ; j++)
            {
                std::cout << j << "\t" << data[i*num_diag+j] << std::endl;
            }
            std::cin.get();
        }*/
        if(nd!=num_diag)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The number of diagonals upon checking was not consisten");
        }
    }

}
