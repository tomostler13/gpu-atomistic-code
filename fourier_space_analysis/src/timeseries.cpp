#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
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
    std::string dump;
    std::string idump;
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
    }
    Array<fftw_complex> Sq;
    //need to read the file and set the variables here




    Sq.resize(num_samples);
    FIXOUT(Info,"Planning the time-series analysis:" << std::flush);
    //At the end we want to perform the FFT in time
    fftw_plan fttime;
    //We're going to use the same identical plan quite a few times so it's
    //worth plannig it well
    fftw_set_timelimit(60);
    fttime=fftw_plan_dft_1d(num_samples,Sq.ptr(),Sq.ptr(),FFTW_FORWARD,FFTW_PATIENT);
    /*SUCCESS(config::Info);
    FIXOUT(config::Info,"Executing the time-series:" << std::flush);
    fftw_execute(fttime);
    SUCCESS(config::Info);
            int negq=-static_cast<int>(i)+num_samples;
            double freq=(static_cast<double>(i)*2.0*M_PI)/(static_cast<double>(dsf::dsfupdate*spins::update*num_samples)*llg::dt);
//calculate the bits of the one-sided PSD
            const double Hq = Cqt(k,i)[0]*Cqt(k,i)[0] + Cqt(k,i)[1]*Cqt(k,i)[1];
            const double Hmq = Cqt(k,negq)[0]*Cqt(k,negq)[0] + Cqt(k,negq)[1]*Cqt(k,negq)[1];
        */
    return(0);
}

