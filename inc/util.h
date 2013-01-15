// File: util.h
// Author:Tom Ostler
// Last-modified: 03 Jan 2013 15:41:30
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <fftw3.h>
#include "../inc/array.h"
#include "../inc/array3d.h"
#ifndef _UTIL_H_
#define _UTIL_H_
namespace util
{
    extern bool ui;
    extern bool pvout;
    extern std::string pvf;
    extern std::string dir;
    //convergence class
	class RunningStat
    {
    public:
        RunningStat() : m_n(0) {}

        void Clear()
        {
            m_n = 0;
        }

        void Push(double x)
        {
            m_n++;

            // See Knuth TAOCP vol 2, 3rd edition, page 232
            if (m_n == 1)
            {
                m_oldM = m_newM = x;
                m_oldS = 0.0;
            }
            else
            {
                m_newM = m_oldM + (x - m_oldM)/m_n;
                m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

                // set up for next iteration
                m_oldM = m_newM;
                m_oldS = m_newS;
            }
        }

        int NumDataValues() const
        {
            return m_n;
        }

        double Mean() const
        {
            return (m_n > 0) ? m_newM : 0.0;
        }

        double PrevMean() const
        {
            return (m_n > 0) ? m_oldM : 0.0;
        }

        double Variance() const
        {
            return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
        }

        double StandardDeviation() const
        {
            return sqrt( Variance() );
        }
        bool is_conv() const
        {
            return(isconv);
        }

    private:
        int m_n;
        bool isconv;
        double m_oldM, m_newM, m_oldS, m_newS;
    };

    extern unsigned int update;
    void initUtil(int argc,char *argv[]);
    void outputSpinsVTU(unsigned int);
    std::string getdir(unsigned int);
    void calcm();
    void fillfloat(int,double*,float*);
    void fillfloat(int,int,int,Array3D<fftw_complex>,Array3D<fftwf_complex>);
    void copy3vecto1(int,double*,double*,double*,double*);
    void copy3vecto1(int,float*,float*,float*,float*);
    void copy3vecto1(int,Array<double>,Array<double>,Array<double>,double*);
}
#endif /*_UTIL_H_*/
