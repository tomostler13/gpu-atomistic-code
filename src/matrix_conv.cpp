//File matrix_conv.cpp
// Author: Tom Ostler
// Created: 07 Jan 2014
// Last-modified: 07 Oct 2014 11:19:28
// The routines within this file convert a 2D matrix to
// a number of formats depending on the routine used. The
// return structures depend on the storage format
#include "../inc/arrays.h"
#include "../inc/matrix_conv.h"
#include <string>
namespace matconv
{

    // This calculate the diagonal offsets for the DIA format
    // The first arguement is the array to store the offsets
    // The second is the J matrix (4 dimensions because it
    // represents all elements of the interactio matrix).
    // This function is overloaded. If the arguement list
    // contains only one offset array then we calculate that
    // for the diagonals of the exchange tensor.
    void dia_offsets(Array<int>& offset,Array4D<double>& JMat,unsigned int& num_diag,unsigned int N)
    {
        //for counting the number of non-zero diagonals
        unsigned int dia_count=0;
        //set the number of diagonals to zero
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

        }
    }
}
