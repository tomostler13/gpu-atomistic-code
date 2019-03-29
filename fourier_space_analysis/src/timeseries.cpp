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
    if(argc<5)
    {
        std::cerr << "Error: You must give a file, the width in Hz as paramters and a min and max frequency for finding the max amplitude." << std::endl;
        exit(0);
    }
    std::string infofile=argv[1];
    double width=atof(argv[2]);
    double fcutl=atof(argv[3]);
    double fcutu=atof(argv[4]);
    std::ofstream Info;
    Info.open("info.dat");
    if(!Info.is_open())
    {
        std::cerr << "Error: Could not open the information (output) file." << std::endl;
        exit(0);
    }
    if(argc < 6)
    {
        FIXOUT(Info,"You have chosen no time window" << std::endl);
    }
    else
    {
        FIXOUT(Info,"You have chosen the window:" << window << std::endl);
    }

    FIXOUT(Info,"Lower frequency cut off:" << fcutl << " [rad/s]" << std::endl);
    FIXOUT(Info,"Upper frequency cut off:" << fcutu << " [rad/s]" << std::endl);
    if(fcutl>fcutu)
    {
        std::cerr << "Lower frequency cut-off cannot be larger than the upper cut-off" << std::endl;
        exit(0);
    }
    if(argc > 6)
    {

        window=argv[5];
        if(window=="hamm")
        {
            windowi=1;
            if(argc > 6 && argc < 8)
            {
                std::cerr << "You have specified a windowing function the use of a generalized hamming window without providing alpha and beta" << std::endl;
                std::cerr << "Usage: ./timeseries <file> <window> <params>" << std::endl;
                std::cerr << "For a generalized Hamming window use, <window>=hamm, params = alpha beta." << std::endl;
            }
            else
            {
                hammalpha=atof(argv[5]);
                hammbeta=atof(argv[6]);
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

    const double normsize=static_cast<double>(kdim[0])*static_cast<double>(kdim[1])*static_cast<double>(kdim[2])*static_cast<double>(num_samples);
    FIXOUT(Info,"Normalization constant" << normsize << std::endl);
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
    fttime=fftw_plan_dft_1d(num_samples,Sq.ptr(),Sq.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
    SUCCESS(Info);

    if(windowi>1)
    {
        std::cerr << "Error: you are trying to use a window that has not been implemented yet" << std::endl;
    }

    //need to convert the width to number of pixels (samples)
    //one pixel corresponds to 2*pi/(num_samples*dt)
    //so one pixel in Hz is 1/(num_samples*dt)
    unsigned int pixelwidth=static_cast<unsigned int>(width*static_cast<double>(num_samples)*dt+0.5);

    FIXOUT(Info,"The desired width of the gaussian blur is:" << width << " [Hz]" << std::endl);
    FIXOUT(Info,"The number of pixels corresponding to our width is:" << pixelwidth << std::endl);
    if(pixelwidth<1)
    {
        std::cerr << "Error: the pixel width for the Gaussian smoothing cannot be less than 1." << std::endl;
        exit(0);
    }
    FIXOUT(Info,"The actual smoothing width (due to rounding) is:" << static_cast<double>(pixelwidth)*static_cast<double>(num_samples)*dt << std::endl);
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
        FIXOUT(Info,"Normalizing data:" << std::flush);
        for(unsigned int t = 0 ; t < num_samples ; t++)
        {
            Sq(t)[0]/=normsize;
            Sq(t)[1]/=normsize;
        }
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
        Array<double> res,ores;
        res.resize(num_samples/2);
        ores.resize(num_samples/2);
        res.IFill(0);
        ores.IFill(0);
        for(unsigned int w = 0 ; w < num_samples/2 ; w++)
        {
        //    std::cout << "w = " << w << std::endl;
            int negq=-static_cast<int>(w)+num_samples-1;
            //calculate the bits of the one-sided PSD
            const double Hq = Sq(w)[0]*Sq(w)[0] + Sq(w)[1]*Sq(w)[1];
        //    std::cout << negq << std::endl;
            const double Hmq = Sq(negq)[0]*Sq(negq)[0] + Sq(negq)[1]*Sq(negq)[1];
            res(w)=Hq+Hmq;
            ores(w)=res(w);
        }
        FIXOUT(Info,"Applying 1D Gaussian filter:" << std::flush);
        gaussianiir1d(ores,res,num_samples/2,pixelwidth,10);
        SUCCESS(Info);
        FIXOUT(Info,"Finding max value of both smoothed and raw data:" << std::endl);
        double maxraw=0.0;
        double maxsmooth=0.0;
        unsigned int maxrawindex=0;
        unsigned int maxsmoothindex=0;
        for(unsigned int w = 0 ; w < num_samples/2 ; w++)
        {
            if(res(w)>maxsmooth && static_cast<double>(w)*2*M_PI/(static_cast<double>(num_samples)*dt) > fcutl && static_cast<double>(w)*2*M_PI/(static_cast<double>(num_samples)*dt) < fcutu)
            {
                maxsmooth=res(w);
                maxsmoothindex=w;
            }
            if(ores(w)>maxraw && w>1 && static_cast<double>(w)*2*M_PI/(static_cast<double>(num_samples)*dt) > fcutl && static_cast<double>(w)*2*M_PI/(static_cast<double>(num_samples)*dt) < fcutu)
            {
                maxraw=ores(w);
                maxrawindex=w;
            }
        }
        SUCCESS(Info);
        FIXOUT(Info,"The maximum value of the raw data is:" << maxraw << std::endl);
        FIXOUT(Info,"This occurs at frequency: " << (static_cast<double>(maxrawindex)*2.0*M_PI)/(static_cast<double>(num_samples)*dt) << " [rad/s]"  << std::endl);
        FIXOUT(Info,"The maximum value of the smoothed data is:" << maxsmooth << std::endl);
        FIXOUT(Info,"This occurs at frequency: " << (static_cast<double>(maxsmoothindex)*2.0*M_PI)/(static_cast<double>(num_samples)*dt) << " [rad/s]" << std::endl);
        for(unsigned int w = 0 ; w < num_samples/2 ; w++)
        {
            double freq=(static_cast<double>(w)*2.0*M_PI)/(static_cast<double>(num_samples)*dt);
//calculate the bits of the one-sided PSD
            dataout << freq << "\t" << ores(w) << "\t" << res(w) << "\t" << ores(w)/maxraw << "\t" << res(w)/maxsmooth << std::endl;
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
