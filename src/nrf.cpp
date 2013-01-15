// File: nrf.cpp
// Author:Tom Ostler
// Last-modified: 18 Dec 2012 17:39:58
//
// Functions from Numerical recipes. Slightly modified in
// places to be C++ compatible.
#include <math.h>
#include <iostream>
#include <stdio.h>
#include "../inc/nrutil.h"
#include "../inc/nrf.h"
#include "../inc/mf.h"
#include "../inc/tdp.h"
#define NRANSI
#define FREERETURN {free_dmatrix(fjac,1,n,1,n);free_dvector(fvec,1,n);\
	free_dvector(p,1,n);free_ivector(indx,1,n);return;}

namespace nrf
{
    void ludcmp(double **a, int n, int *indx, double& d)
    {
        int imax;
        double big=0.0,dum=0.0,sum=0.0,temp=0.0;
        double *vv;
        //vv stores the implicit scaling of each row
        vv=dvector(1,n);//this declaration requires the NR utilities
        d=1.0;//No row interchanges yet.
        for (int i=1;i<=n;i++)
        {
            //Loop over rows to get the implicit scaling information
            big=0.0;
            for(int j=1;j<=n;j++)
            {
                if ((temp=fabs(a[i][j])) > big) big=temp;
            }
            //No nonzero largest element
            if (big == 0.0)
            {
                std::cerr << "Error: singular matrix in routine ludcmp" << std::endl;
            }
            vv[i]=1.0/big;//save the scaling
        }//end of i loop
        for (int j=1;j<=n;j++)
        {
            //This is the loop over the columns of Crout's method.
            for (int i=1;i<j;i++)
            {
                //This is equation (2.3.12) except for i=j from NR in C
        //printf("\ni=%d",i);
                sum=a[i][j];
                for (int k=1;k<i;k++)
                {     /*printf("\nk=%d",k);*/
                    sum -= a[i][k]*a[k][j];
                }
                    a[i][j]=sum;
            }
            big=0.0;//Initialize for the search for largest pivot element.
            for (int i=j;i<=n;i++)
            {
                //This is the i=j part of equation (2.3.12) and i = k +1 . . . N, of equation (2.3.12).
                sum=a[i][j];
                for (int k=1;k<j;k++)
                {
                    sum -= a[i][k]*a[k][j];
                }
                a[i][j]=sum;
                //Is the figure of merit for the pivot better than the best so far?
                if ( (dum=vv[i]*fabs(sum)) >= big)
                {
                    big=dum;
                    imax=i;
                }
            }
            if (j != imax)//Do we need to interchange rows?
            {
                //printf("%d\t%d\n",j,imax);
                for (int k=1;k<=n;k++)
                {
                    //Yes, do so...
                    dum=a[imax][k];
                    a[imax][k]=a[j][k];
                    a[j][k]=dum;
                }
                d = -(d);//...and change the parity of d.
                vv[imax]=vv[j];//Also interchange the scale factor.
            }
            indx[j]=imax;
            if (a[j][j] == 0.0)
            {
                a[j][j]=1e-40;
            }
            //If the pivot element is zero the matrix is singular (at least to the precision of the
            //algorithm). For some applications on singular matrices, it is desirable to substitute
            //TINY for zero.

            if (j != n)
            {
            //printf("Strangely satisfied...%d\t%d\n",j,n);
                dum=1.0/(a[j][j]);
                for (int i=j+1;i<=n;i++)
                {
                    a[i][j] *= dum;
                }
            }
        }//go back for the next column in the reduction
        free_dvector(vv,1,n);
    }
    //the arguements are as follows:
    //ntrail - number of trial steps
    //x - the sublattice magnetizations
    //n - the number of sublattices
    //tolx - convergence parameter
    //tolf - function sum convergence parameter
    void mnewt(int& ntrial, double *x,int n, double& tolx, double& tolf,unsigned int& atom)
    {
        mf::beta[atom] = 1./(1.38e-16*tdp::systemp[atom]);
        int k,i,*indx;
        double errx,errf,d,*fvec,**fjac,*p;

        indx=ivector(1,n);
        for(unsigned int i = 1 ; i <= n ; i++)
        {
            indx[i]=0;
        }
        p=dvector(1,n);
        fvec=dvector(1,n);
        fjac=dmatrix(1,n,1,n);
        for (k=1;k<=ntrial;k++) {
            mf::usrfun(x,n,fvec,fjac,mf::beta[atom],atom);
            errf=0.0;
            for (i=1;i<=n;i++) errf += fabs(fvec[i]);
            if (errf <= tolf) FREERETURN
            for (i=1;i<=n;i++) p[i] = -fvec[i];
            ludcmp(fjac,n,indx,d);
            lubksb(fjac,n,indx,p);
            errx=0.0;
            for (i=1;i<=n;i++) {
                errx += fabs(p[i]);
                x[i] += p[i];
            }
            if (errx <= tolx) FREERETURN
        }
        FREERETURN

        #undef NRANSI
        #undef FREERETURN
    }
    void lubksb(double **a, int n, int *indx, double b[])
    {
        int i,ii=0,ip,j;
        double sum=0.0;

        for (i=1;i<=n;i++) {
            ip=indx[i];
            sum=b[ip];
            b[ip]=b[i];
            if (ii)
                for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
            else if (sum) ii=i;
            b[i]=sum;
        }
        for (i=n;i>=1;i--) {
            sum=b[i];
            for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
            b[i]=sum/a[i][i];
        }
    }
    /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

}
