// File: nrf.h
// Author:Tom Ostler
// Last-modified: 17 Dec 2012 10:30:15
//
// This header file contains extern declarations and
// function prototypes for any Numerical recipes functions
// used in this project.
//
#include "../inc/array.h"
#ifndef _NRF_H_
#define _NRF_H_
namespace nrf
{
    void ludcmp(double **a,int n,int *indx,double& d);
    void mnewt(int& ntrial,double *x,int n,double& tolx,double& tolf,unsigned int& atom);
    void lubksb(double **a,int n,int *indx,double b[]);

}
#endif /*_NRF_H_*/
