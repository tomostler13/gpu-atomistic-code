#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include "../../inc/array.h"
#include "../../inc/array2d.h"
#include "../../inc/defines.h"
//function prototype for Gaussian smoothing
void gaussianiir1d(Array<double>&,Array<double>&, unsigned int, double,unsigned int);
//This code is designed to work with the atomistic spin
//dynamics code and post-process the fourier components S(\vec{q},t)
//for a number of the vectors of q. The atomistic spin dynamics code
//provides two files 1) and information file that contains the info
//on the second file, 2) the data of the fourier transformed spin
//map for given q-values. The 3D FFT of the spin map has already
//been performed during the integration of the LLG equation. To
//perform a time series (i.e. to calculate the DSF) you should use
//the timeseries code. Otherwise, what this code does is multiply
//the real and imaginary components of the fourier transformed spin
//map (for each fourier component) and put the information together
//as time-sliced data. The data is output as both raw, normalized
//raw, smoothed and smooted normalized
int main(int argc,char *argv[])
{
    double idump;
    std::string dump;
    if(argc<4)
    {
        std::cerr << "Error: You must give a file and the width of the smoothing in reciprocal lattice spacings as parameters and the number of lines\n./ssf <file> <width> <num_timesteps>" << std::endl;
    }
    std::string infofile=argv[1];
    double width=atof(argv[2]);
    unsigned int pixelwidth=width;

    unsigned int num_ts=atoi(argv[3]);
    std::ofstream Info;
    Info.open("info.dat");
    if(!Info.is_open())
    {
        std::cerr << "Error: Could not open the information (output) file." << std::endl;
        exit(0);
    }
    FIXOUT(Info,"Number of timesteps:" << num_ts << std::endl);
    FIXOUT(Info,"Width of smoothing function:" << width << " [reciprocal lattice spacings]" << std::endl);
    FIXOUT(Info,"K-vec information file:" << infofile << std::endl);
    std::ifstream infoin;
    infoin.open(infofile.c_str());
    if(infoin.is_open())
    {
        FIXOUT(Info,"K-vec information file opened" << std::endl);
    }
    else
    {
        std::cerr << "Error: Could not open the k-vector information file." << std::endl;
    }

    std::string datfile;
    infoin >> dump;
    infoin >> datfile;
    FIXOUT(Info,"Data file:" << datfile << std::endl);
    int num_samples=0,dim[3]={0,0,0},kdim[3]={0,0,0},nk=0;
    double dt=0.0;
    infoin >> dump;
    infoin >> dim[0] >> dim[1] >> dim[2];
    FIXOUTVEC(Info,"Dimensions of original system (unit cells):",dim[0],dim[1],dim[2]);
    infoin >> dump;
    infoin >> kdim[0] >> kdim[1] >> kdim[2];
    FIXOUTVEC(Info,"No. K-points in each direction:",kdim[0],kdim[1],kdim[2]);
    infoin >> dump;
    infoin >> nk;
    unsigned int maxk=nk;
    infoin >> dump;
    unsigned ignorenum=0;
    Array<int> ignorelist;
    bool ignoresome=false;
    if(dump=="#Ignore")
    {
        //then we have to read a set of k-point to ignore
        infoin >> ignorenum;
        //decrease the number of k-points
        nk-=ignorenum;
        ignoresome=true;
        ignorelist.resize(ignorenum);
        ignorelist.IFill(0);
        FIXOUT(Info,"We are ignorging some k-vectors for output (how many):" << ignorenum << std::endl);
        FIXOUT(Info,"List of ignored vectors:" << "[ " << std::flush);
        for(unsigned int i = 0 ; i < ignorenum ; i++)
        {
            infoin >> ignorelist[i];
            if(i<(ignorenum-1))
            {
                Info << ignorelist[i] << " , ";
            }
            else
            {
                Info << ignorelist[i] << " ]" << std::endl;
            }
        }
    }
    else
    {
        FIXOUT(Info,"Not ignoring any k-vectors" << std::endl);
    }

    FIXOUT(Info,"Number of k points sampled:" << nk << std::endl);
    const double normsize=static_cast<double>(kdim[0])*static_cast<double>(kdim[1])*static_cast<double>(kdim[2]);
    FIXOUT(Info,"Normalization constant (per timesteps):" << normsize << std::endl);
    Array2D<int> kvecs;
    kvecs.resize(nk,3);
    unsigned int iglu=0;
    for(unsigned int i = 0 ; i < maxk ; i++)
    {
        if(i==0 && ignoresome==false)
        {
            infoin >> idump >> dump;
        }
        else if(i>0 && ignoresome==false)
        {
            infoin >> dump >> idump >> dump;
        }
        else
        {
             infoin >> dump >> idump >> dump;
        }
        int tk[3]={0,0,0};
        infoin >> tk[0] >> tk[1] >> tk[2];//kvecs(i,0) >> kvecs(i,1) >> kvecs(i,2);
        if(ignoresome==true)
        {
            if(i==ignorelist[iglu])
            {
                iglu++;
            }
            else
            {
                std::stringstream sstr;
                sstr << "K-vector " << i << " of " << nk << ":";
                std::string str=sstr.str();
                kvecs(i,0)=tk[0];
                kvecs(i,1)=tk[1];
                kvecs(i,2)=tk[2];
                FIXOUTVEC(Info,str,kvecs(i,0),kvecs(i,1),kvecs(i,2));

            }
        }
        else
        {
            std::stringstream sstr;
            sstr << "K-vector " << i << " of " << nk << ":";
            std::string str=sstr.str();
            kvecs(i,0)=tk[0];
            kvecs(i,1)=tk[1];
            kvecs(i,2)=tk[2];
            FIXOUTVEC(Info,str,kvecs(i,0),kvecs(i,1),kvecs(i,2));
        }
        infoin >> dump >> idump >> dump >> idump >> idump >> idump;
    }
    //allocate an array to store the S(q) at one instance in time
    Array<double> S,oS;
    S.resize(nk);oS.resize(nk);
    S.IFill(0);oS.IFill(0);
    std::ifstream datain(datfile.c_str());
    if(!datain.is_open())
    {
        std::cerr << "Error: Could not open file for reading data"<< std::endl;
        exit(0);
    }
    std::ofstream dataout("ssf.dat");
    if(!dataout.is_open())
    {
        std::cerr << "Error: Could not open data output file" << std::endl;
        exit(0);
    }
    for(unsigned int t = 0 ; t < num_ts ; t++)
    {
        S.IFill(0);
        oS.IFill(0);
        //real and imaginary parts
        double Sr=0.0,Sc=0.0;
        //the time index
        double timein=0.0;
        datain >> timein;
        unsigned int ig=0;
        unsigned int kcount=0;
        double maxraw=0.0;
        unsigned int maxrawindex=0;
        double maxsmooth=0.0;
        unsigned int maxsmoothindex=0;
        //read the file
        for(unsigned int i = 0 ; i < maxk ; i++)
        {
            if(ignoresome==true)
            {
                if(i==ignorelist[ig])
                {
                    datain >> Sr;
                    datain >> Sc;
                }
                else
                {
                    datain >> Sr;
                    datain >> Sc;
                    S[kcount]=Sr*Sr+Sc*Sc;
                    S[kcount]/=normsize;
                    if(S[kcount]>maxraw)
                    {
                        maxraw=S[kcount];
                        maxrawindex=kcount;
                    }
                    kcount++;
                }
            }
            else
            {
                datain >> Sr;
                datain >> Sc;
                S[kcount]=Sr*Sr+Sc*Sc;
                S[kcount]/=normsize;
                if(S[kcount]>maxraw)
                {
                    maxraw=S[kcount];
                    maxrawindex=kcount;
                }
                kcount++;
            }
        }

        //smooth the data
        gaussianiir1d(S,oS,nk,pixelwidth,10);
        //find max of smoothed data
        for(unsigned int i = 0 ; i < nk ; i++)
        {
            if(oS[i]>maxsmooth)
            {
                maxsmooth=oS[i];
            }
        }
        for(unsigned int i = 0 ; i < nk ; i++)
        {

            dataout << timein << "\t" << kvecs(i,0) << "\t" << kvecs(i,1) << "\t" << kvecs(i,2) << "\t" << S[i] << "\t" << S[i]/maxraw << "\t" << oS[i] << "\t" << oS[i]/maxsmooth << std::endl;
        }
        dataout << std::endl << std::endl;
    }
    dataout.close();
    if(dataout.is_open())
    {
        std::cerr << "Could not close the output data file." << std::endl;
    }
    datain.close();
    if(datain.is_open())
    {
        std::cerr << "Could not close the input data file." << std::endl;
    }

    return(0);
}

/*
 * Reference:
 * Alvarez, Mazorra, "Signal and Image Restoration using Shock Filters and
 * Anisotropic Diffusion," SIAM J. on Numerical Analysis, vol. 31, no. 2,
 * pp. 590-605, 1994.
*/
void gaussianiir1d(Array<double>& indata,Array<double>& data, unsigned int length, double sigma,unsigned int numsteps)
{
    double lambda, dnu;
    double nu, boundaryscale, postscale;
    long i;
    int step;
    if(!data.ptr() || length < 1 || sigma <= 0)
        return;
    lambda = (sigma*sigma)/(2.0*numsteps);
    dnu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
    nu = (double)dnu;
    boundaryscale = (double)(1.0/(1.0 - dnu));
    postscale = (double)(pow(dnu/lambda,numsteps));
    //copy to indata to preserve the original data (we are not finished with it)
    for(unsigned int i = 0 ; i < length ; i++)
    {
        data[i]=indata[i];
    }
    for(step = 0; step < numsteps; step++)
    {
        data[0] *= boundaryscale;

        /* Filter rightwards (causal) */
        for(i = 1; i < length; i++)
            data[i] += nu * data[i - 1];

        data[i = length - 1] *= boundaryscale;

        /* Filter leftwards (anti-causal) */
        for(; i > 0; i--)
            data[i - 1] += nu*data[i];
    }

    for(i = 0; i < length; i++)
        data[i] *= postscale;

    return;
}
