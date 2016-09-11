#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "array.h"
#include "array2d.h"
#include "array3d.h"
#include "array4d.h"
int main()
{
    //Number of neighbours
    int nJij=6;
    //Number of species
    const int nspec=2;
    //Magnetic Moments
    Array<double> mu;
    mu.resize(nspec);
    mu.IFill(0.0);
    mu[0]=2.00;
    mu[1]=2.00;
    //spin directions
    const double s[2]={0.999,0.999};
    //amounts of each sublattice
    const double x[2]={0.9999999,0.0000001};
    //some constants
    const double muB=9.27e-24;
    const double gamma=1.76e11;
    double a=300e-12;
    //for outputting the elements of various stuff
    std::ofstream *ofs=NULL;ofs=new std::ofstream[nspec];
    std::ofstream test("lswt.dat");
    for(unsigned int i = 0 ; i < nspec ; i++)
    {
        std::stringstream sstr;sstr << "LSWT_spec_" << i << ".dat";
        std::string str=sstr.str();
        ofs[i].open(str.c_str());
        if(!ofs[i].is_open())
        {
            std::cerr << "Could not open file " << str << std::endl;
        }
    }
    std::ifstream *inrJ=NULL;inrJ = new std::ifstream[nspec*nspec];
    for(unsigned int i = 0 ; i < nspec ; i++)
    {
        for(unsigned int j = 0 ; j < nspec ; j++)
        {
            std::stringstream sstr;sstr << "rJ" << i << j << ".dat";
            std::string str=sstr.str();
            inrJ[i*nspec+j].open(str.c_str());
            if(!inrJ[i*nspec+j].is_open())
            {
                std::cerr << "Error opening file " << str << " for reading." << std::endl;
                exit(0);
            }
        }
    }
    Array4D<double> rvec;
    rvec.resize(nJij,nspec,nspec,3);
    rvec.IFill(0);
    Array3D<double> Jij;
    Jij.resize(nJij,nspec,nspec);
    Jij.IFill(0);
    double dump;
    //read in the exchange
    for(unsigned int i = 0 ; i < nJij ; i++)
    {
        for(unsigned int s1 = 0 ; s1 < nspec ; s1++)
        {
            for(unsigned int s2 = 0 ; s2 < nspec ; s2++)
            {
                inrJ[s1*nspec+s2] >> rvec(i,s1,s2,0) >> rvec(i,s1,s2,1) >> rvec(i,s1,s2,2) >> Jij(i,s1,s2);
                //convert from eV to J
                Jij(i,s1,s2)*=1.6e-19;
            }
        }
    }
    for(unsigned int i = 0 ; i < nspec ; i++)
    {
        for(unsigned int j = 0 ; j < nspec ; j++)
        {
            inrJ[i*nspec+j].close();
            if(inrJ[i*nspec+j].is_open())
            {
                std::cerr << "Error closing input file between species " << i << " and " << j << std::endl;
            }
        }
    }
    //calculate the dispersion for the individual species
    //Gamma to X
    Array<double> Rsum,Csum;
    Rsum.resize(nspec);
    Csum.resize(nspec);
    Rsum.IFill(0);
    Csum.IFill(0);
    //In this case kx=0,kz=0
    double k[3]={0,0,0};
    for(k[1] = 0.0 ; k[1] < 2*M_PI/a ; k[1]+=((2.0*M_PI/a)/500.0))
    {
        double J00sum=0.0,J11sum=0.0,J01sum=0.0,J10sum=0.0;
        double M[nspec][nspec] = {{0.0,0.0},{0.0,0.0}};
        //loop over all of the R vectors
        for(unsigned int r = 0 ; r < nJij ; r++)
        {
            //This loop is ONLY for the pure interactions (i.e. if the two species were bulk)
            for(unsigned int i = 0 ; i < nspec ; i++)
            {
                const double R[3]={rvec(r,i,i,0)*0.5*a,rvec(r,i,i,1)*0.5*a,rvec(r,i,i,2)*0.5*a};
                const double dotp=k[0]*R[0]+k[1]*R[1]+k[2]*R[2];
                //for each R calculate the frequency
                Rsum[i]+=((gamma/(mu[i]*9.27e-24))*Jij(r,i,i)*(1-cos(dotp)));
                Csum[i]+=((gamma/(mu[i]*9.27e-24))*Jij(r,i,i)*(sin(dotp)));
            }
            const double R00[3]={rvec(r,0,0,0)*0.5*a,rvec(r,0,0,1)*0.5*a,rvec(r,0,0,2)*0.5*a};
            const double R01[3]={rvec(r,0,1,0)*0.5*a,rvec(r,0,1,1)*0.5*a,rvec(r,0,1,2)*0.5*a};
            const double R10[3]={rvec(r,1,0,0)*0.5*a,rvec(r,1,0,1)*0.5*a,rvec(r,1,0,2)*0.5*a};
            const double R11[3]={rvec(r,1,1,0)*0.5*a,rvec(r,1,1,1)*0.5*a,rvec(r,1,1,2)*0.5*a};
            const double dotp00=k[0]*R00[0]+k[1]*R00[1]+k[2]*R00[2];
            const double dotp01=k[0]*R01[0]+k[1]*R01[1]+k[2]*R01[2];
            const double dotp10=k[0]*R10[0]+k[1]*R10[1]+k[2]*R10[2];
            const double dotp11=k[0]*R11[0]+k[1]*R11[1]+k[2]*R11[2];
            //1-x -> x[0]
            //x -> x[1]
            M[0][0] += (gamma/(mu[0]*9.27e-24))*x[0]*Jij(r,0,0)*(1.0-cos(dotp00))*s[0] + 0.5*(gamma/(mu[0]*9.27e-24))*x[1]*Jij(r,0,1)*s[1] + 0.5*(gamma/(mu[1]*9.27e-24))*x[0]*Jij(r,1,0)*s[0];
            M[0][1] += (gamma/(mu[0]*9.27e-24))*x[1]*Jij(r,0,1)*(cos(dotp01))*s[1];
            M[1][0] += (gamma/(mu[1]*9.27e-24))*x[0]*Jij(r,1,0)*(cos(dotp10))*s[0];
            M[1][1] += (gamma/(mu[1]*9.27e-24))*x[1]*Jij(r,1,1)*(1.0-cos(dotp11))*s[1] + 0.5*(gamma/(mu[1]*9.27e-24))*x[0]*Jij(r,1,0)*s[0] + 0.5*(gamma/(mu[0]*9.27e-24))*x[1]*Jij(r,0,1)*s[1];
            J00sum+=Jij(r,0,0);
            J11sum+=Jij(r,1,1);
            J01sum+=Jij(r,0,1);
            J10sum+=Jij(r,1,0);

        }

        for(unsigned int i = 0 ; i < nspec ; i++)
        {
            ofs[i] << k[1]*a << "\t" << Rsum[i] << "\t" << Csum[i] << std::endl;
        }
        test << k[1]*a << "\t" << M[0][0] << "\t" << M[0][1] << "\t" << M[1][0] << "\t" << M[1][1] << "\t" << 0.5*(sqrt( (M[0][0]+M[1][1])*(M[0][0]+M[1][1]) - 4.0*(M[0][0]*M[1][1]-M[0][1]*M[1][0]))+M[1][1]+M[0][0]) << "\t" << 0.5*(-sqrt( (M[0][0]+M[1][1])*(M[0][0]+M[1][1]) - 4.0*(M[0][0]*M[1][1]-M[0][1]*M[1][0]))+M[0][0]+M[1][1]) << "\t" << (0.5*( ((M[0][0]-M[1][1])/(sqrt((M[0][0]+M[1][1])*(M[0][0]+M[1][1])-4.0*M[0][1]*M[1][0])))+1   )) << "\t" << (0.5*( ((M[0][0]-M[1][1])/(sqrt((M[0][0]+M[1][1])*(M[0][0]+M[1][1])-4.0*M[0][1]*M[1][0]))) -1  )) << std::endl;
        for(unsigned int i = 0 ; i < nspec ; i++)
        {
            Rsum[i]=0.0;
            Csum[i]=0.0;
        }
        std::cout << J00sum << "\t" << J01sum << "\t" << J10sum << "\t" << J11sum << std::endl;
        //exit(0);
    }
    for(unsigned int i = 0 ; i < nspec ; i++)
    {
        ofs[i].close();
        if(ofs[i].is_open())
        {
            std::cerr << "Error closing output file, id=" << i << std::endl;
        }
    }
    return(0);
}
