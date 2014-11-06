#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include "../../inc/array.h"
#include "../../inc/array2d.h"
#include "../../inc/defines.h"
//This code is designed to work with the atomistic spin
//dynamics code and post-process the fourier components S(\vec{q},t)
//for a number of the vectors of q. The atomistic spin dynamics code
//provides two files 1) and information file that contains the info
//on the second file, 2) the data of the fourier transformed spin
//map for given q-values. Whilst fftw_plan_many provides a more elegant
//means of performing the analysis the amount of memory required can
//grow quite a lot if many q-vectors with a long time-series are
//analysed. To that end the file (2) is re-read for each q-vector. This
//could be quite time consuming but overcomes the memory issue.
//The second filename is taken from the info file (1).
int main(int argc,char *argv[])
{
    double idump;
    double hammalpha=0.0,hammbeta=0.0;
    std::string dump;
    std::string window;
    //this tells us which window we are using
    //0 - no window
    //1 - generalized hamming
    unsigned int windowi=0;
    if(argc<2)
    {
        std::cerr << "Error: You must give a file as an arguement." << std::endl;
        exit(0);
    }
    std::string infofile=argv[1];
    std::ofstream Info;
    Info.open("info.dat");
    if(!Info.is_open())
    {
        std::cerr << "Error: Could not open the information (output) file." << std::endl;
        exit(0);
    }
    if(argc < 3)
    {
        FIXOUT(Info,"You have chosen no time window" << std::endl);
    }
    else
    {
        FIXOUT(Info,"You have chosen the window:" << window << std::endl);
    }
    if(argc > 3)
    {

        window=argv[2];
        if(window=="hamm")
        {
            windowi=1;
            if(argc > 3 && argc < 5)
            {
                std::cerr << "You have specified a windowing function the use of a generalized hamming window without providing alpha and beta" << std::endl;
                std::cerr << "Usage: ./timeseries <file> <window> <params>" << std::endl;
                std::cerr << "For a generalized Hamming window use, <window>=hamm, params = alpha beta." << std::endl;
            }
            else
            {
                hammalpha=atof(argv[2]);
                hammbeta=atof(argv[3]);
                FIXOUT(Info,"Hamm alpha:" << hammalpha << std::endl);
                FIXOUT(Info,"Hamm beta:" << hammbeta << std::endl);
            }
        }
        else
        {
            std::cerr << "Windowing technique not recognized" << std::endl;
        }
    }
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
    infoin >> num_samples;
    FIXOUT(Info,"Number of samples:" << num_samples << std::endl);
    infoin >> dump;
    infoin >> dt;
    FIXOUT(Info,"Sample timestep:" << dt << std::endl);
    infoin >> dump;
    infoin >> dim[0] >> dim[1] >> dim[2];
    FIXOUTVEC(Info,"Dimensions of original system (unit cells):",dim[0],dim[1],dim[2]);
    infoin >> dump;
    infoin >> kdim[0] >> kdim[1] >> kdim[2];
    FIXOUTVEC(Info,"No. K-points in each direction:",kdim[0],kdim[1],kdim[2]);
    infoin >> dump;
    infoin >> nk;
    FIXOUT(Info,"Number of k points sampled:" << nk << std::endl);
    Array2D<int> kvecs;
    kvecs.resize(nk,3);
    for(unsigned int i = 0 ; i < nk ; i++)
    {
        infoin >> dump >> idump >> dump;
        infoin >> kvecs(i,0) >> kvecs(i,1) >> kvecs(i,2);
        std::stringstream sstr;
        sstr << "K-vector " << i << " of " << nk << ":";
        std::string str=sstr.str();
        FIXOUTVEC(Info,str,kvecs(i,0),kvecs(i,1),kvecs(i,2));
        infoin >> dump >> idump >> dump >> idump >> idump >> idump;
    }
    Array<fftw_complex> Sq;
    Sq.resize(num_samples);
    FIXOUT(Info,"Planning the time-series analysis:" << std::flush);
    //At the end we want to perform the FFT in time
    fftw_plan fttime;
    //We're going to use the same identical plan quite a few times so it's
    //worth plannig it well
    fftw_set_timelimit(60);
    fttime=fftw_plan_dft_1d(num_samples,Sq.ptr(),Sq.ptr(),FFTW_FORWARD,FFTW_PATIENT);
    SUCCESS(Info);

    if(windowi>1)
    {
        std::cerr << "Error: you are trying to use a window that has not been implemented yet" << std::endl;
    }

    //loop over the number of k-vectors that we want to perform the time series for
    for(unsigned int i = 0 ; i < nk ; i++)
    {
        std::ifstream datain(datfile.c_str());
        if(!datain.is_open())
        {
            std::cerr << "Error: Could not open file for reading k-vector " << i << std::endl;
            exit(0);
        }
        Sq.IFill(0);
        for(unsigned int t = 0 ; t < num_samples ; t++)
        {
            //get rid of the time index
            datain >> dump;
            //loop over the colums and get rid of the ones we don't want
            for(unsigned int c = 0 ; c < nk ; c++)
            {
                if(c==i)//then we want this column
                {
                    datain >> Sq(t)[0];
                    datain >> Sq(t)[1];
                    if(windowi==1)
                    {
                        Sq(t)[0]*=(hammalpha-hammbeta*cos(2.0*M_PI*static_cast<double>(t)/(static_cast<double>(num_samples)-1.0)));
                        Sq(t)[1]*=(hammalpha-hammbeta*cos(2.0*M_PI*static_cast<double>(t)/(static_cast<double>(num_samples)-1.0)));
                    }
                }
                else
                {
                    datain >> dump >> dump;
                }
            }
        }
        FIXOUTVEC(Info,"Executing fourier transform for k-vector:",kvecs(i,0),kvecs(i,1),kvecs(i,2));
        fftw_execute(fttime);
        SUCCESS(Info);
        //output data
        std::stringstream sstr;
        sstr << "kx" << kvecs(i,0) << "ky" << kvecs(i,1) << "kz" << kvecs(i,2) << ".dat";
        std::string str=sstr.str();
        datain.close();
        if(datain.is_open())
        {
            std::cerr << "Error: Could not close data file on loop " << i << std::endl;
        }

        FIXOUT(Info,"Outputting data to file:" << str << std::endl);
        std::ofstream dataout(str.c_str());
        if(!dataout.is_open())
        {
            std::cerr << "Error: Could not open data file (" << str << ") for writing PSD" << std::endl;
            exit(0);
        }
        for(unsigned int w = 0 ; w < num_samples/2 ; w++)
        {
            int negq=-static_cast<int>(w)+num_samples;
            double freq=(static_cast<double>(w)*2.0*M_PI)/(static_cast<double>(num_samples)*dt);
//calculate the bits of the one-sided PSD
            const double Hq = Sq(w)[0]*Sq(w)[0] + Sq(w)[1]*Sq(w)[1];
            const double Hmq = Sq(negq)[0]*Sq(negq)[0] + Sq(negq)[1]*Sq(negq)[1];
            dataout << freq << "\t" << Hq << "\t" << Hmq << std::endl;
        }

        dataout.close();
        if(dataout.is_open())
        {
            std::cerr << "Error: Could not close the output data file." << std::endl;
            exit(0);
        }
    }
    SUCCESS(Info);
    return(0);
}

