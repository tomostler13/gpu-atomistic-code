// File: util.h
// Author:Tom Ostler
// Last-modified: 17 Jan 2013 14:26:02
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
    void inverse(double*, int);
    void cpuConvFourier();

}
#endif /*_UTIL_H_*/
