// File: util.h
// Author:Tom Ostler
// Last-modified: 09 Dec 2014 20:16:58
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <fftw3.h>
#include <iostream>
#include "../inc/array.h"
#include "../inc/array3d.h"
#ifndef _UTIL_H_
#define _UTIL_H_
namespace util
{
    extern std::ofstream ofs,sofs;
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
    void hcpuConvFourier();
    void dipcpuConvFourier();
    ////////////////////////////////////////////////////////////////////////////////
    //! Compute sum reduction on CPU
    //! We use Kahan summation for an accurate sum of large arrays.
    //! http://en.wikipedia.org/wiki/Kahan_summation_algorithm
    //!
    //! @param data       pointer to input data
    //! @param size       number of input data elements
    ////////////////////////////////////////////////////////////////////////////////
    template<class T>
        T reduceCPU(Array<T> data, unsigned int size)
        {
            T sum = data[0];
            T c = (T)0.0;
            for (int i = 1; i < size; i++)
            {
                T y = data[i] - c;
                T t = sum + y;
                c = (t - sum) - y;
                sum = t;
            }
            return sum;
        }
    double reduceArrayDouble(Array<double>,unsigned int);
    void fillfloat(int,double*,float*);
    void fillfloat(int,int,int,Array3D<fftw_complex>,Array3D<fftwf_complex>);
    void copy3vecto1(int,double*,double*,double*,double*);
    void copy3vecto1(int,float*,float*,float*,float*);
    void copy3vecto1(int,Array<double>,Array<double>,Array<double>,double*);
    std::string exec(char*);
    void outputSpinsVTU(unsigned int);
    void calc_mag();
    void output_mag(unsigned int);
    void init_output();

}
#endif /*_UTIL_H_*/
