#include "array.h"
#include "array2d.h"
#ifndef _MF_H_
#define _MF_H_
namespace mf
{
    extern Array2D<double> J0;
    extern Array<double> fmsub,D;
    extern Array2D<double> I;
    extern double z;
    extern Array2D<double> hmfg;
    extern Array<double> beta,hmf;
    extern bool mfi;

    void initmf(int argc,char *argv[]);

    void usrfun(double *msub,int n,double *fvec,double **fjac,double& beta,unsigned int& atom);
    double langevin(double);
    double dlangevin(double);
}
#endif /*_MF_H_*/
