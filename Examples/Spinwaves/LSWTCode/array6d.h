// File: array5d.h
// Author:Joe Barker
// Last-modified: 20 May 2014 17:00:10
#ifndef __ARRAY6D_H__
#define __ARRAY6D_H__

#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

#include <vector>
#include <cassert>
#include <cstddef>
#include <fftw3.h>

template <typename _Tp>
class Array6D
{
    public:
        typedef unsigned int size_type;
        typedef _Tp value_type;
        typedef _Tp* iterator;
        typedef const _Tp* const_iterator;
        typedef ptrdiff_t difference_type;

        Array6D() : dim0(0), dim1(0), dim2(0), dim3(0), dim4(0), dim5(0), data(0) {}

        Array6D(size_type d0, size_type d1, size_type d2, size_type d3, size_type d4, size_type d5)
            : dim0(d0), dim1(d1), dim2(d2), dim3(d3), dim4(d4), dim5(d5), data(d0*d1*d2*d3*d4*d5) {}

        ~Array6D(){data.clear();}

        inline void clear() {
            dim0 = 0;
            dim1 = 0;
            dim2 = 0;
            dim3 = 0;
            dim4 = 0;
            dim5 = 0;
            data.clear();
        }

        inline int getarrayelement(unsigned int i,unsigned int j,unsigned int k,unsigned int l,unsigned int m,unsigned int n)
        {
            return(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n);
        }

        inline void resize(const size_type d0, const size_type d1, const size_type d2, const size_type d3, const size_type d4,const size_type d5) {
            dim0 = d0; dim1 = d1; dim2 = d2; dim3 = d3; dim4 = d4; dim5=d5;
            data.resize(d0*d1*d2*d3*d4*d5);
        }
        inline void IFill(size_type a)
        {
            for(unsigned int i = 0 ; i < dim0 ; i++)
            {
                for(unsigned int j = 0 ; j < dim1 ; j++)
                {
                    for(unsigned int k = 0 ; k < dim2 ; k++)
                    {
                        for(unsigned int l = 0 ; l < dim3 ; l++)
                        {
                            for(unsigned int m = 0 ; m < dim4 ; m++)
                            {
                                for(unsigned int n = 0 ; n < dim5 ; n++)
                                {
                                    data[((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n]=a;
                                }
                            }
                        }
                    }
                }
            }
        }

        inline _Tp& RESTRICT operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n) {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            return data[((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n];
        }

        inline const _Tp& operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n) const {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            return data[((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n];
        }


        inline _Tp* RESTRICT ptr() {
            return &(data[0]);
        }

    private:
        size_type dim0;
        size_type dim1;
        size_type dim2;
        size_type dim3;
        size_type dim4;
        size_type dim5;
        std::vector<_Tp> data;
};

template <>
class Array6D<fftw_complex>
{
    public:
        typedef unsigned int size_type;
        typedef fftw_complex value_type;
        typedef fftw_complex* iterator;
        typedef const fftw_complex* const_iterator;
        typedef ptrdiff_t difference_type;

        Array6D() : dim0(0), dim1(0), dim2(0), dim3(0), dim4(0), dim5(0), data(0) {}

        Array6D(size_type d0, size_type d1, size_type d2, size_type d3, size_type d4, size_type d5)
        {resize(d0,d1,d2,d3,d4,d5);}

        ~Array6D(){clear();}

        inline void clear() {
            dim0 = 0;
            dim1 = 0;
            dim2 = 0;
            dim3 = 0;
            dim4 = 0;
            dim5 = 0;
            fftw_free(data);
            data = NULL;
        }

        inline void resize(const size_type d0, const size_type d1, const size_type d2, const size_type d3, const size_type d4,const size_type d5) {
            dim0 = d0; dim1 = d1; dim2 = d2; dim3 = d3; dim4 = d4; dim5 = d5;
            if(data != NULL){
                fftw_free(data);
            }
            data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*d0*d1*d2*d3*d4*d5);
        }
        inline void IFill(size_type a)
        {
            for(unsigned int i = 0 ; i < dim0 ; i++)
            {
                for(unsigned int j = 0 ; j < dim1 ; j++)
                {
                    for(unsigned int k = 0 ; k < dim2 ; k++)
                    {
                        for(unsigned int l = 0 ; l < dim3 ; l++)
                        {
                            for(unsigned int m = 0 ; m < dim4 ; m++)
                            {
                                for(unsigned int n = 0 ; n < dim5 ; n++)
                                {
                                    data[((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n][0]=a;
                                    data[((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n][1]=a;
                                }
                            }
                        }
                    }
                }
            }
        }

        inline fftw_complex& RESTRICT operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n) {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            return data[((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n];
        }

        inline const fftw_complex& operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n) const {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            return data[((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n];
        }

        inline fftw_complex* __restrict__ ptr() {
            return data;
        }

        //    inline size_type size(const size_type i) const { return dim[i]; }

    private:
        size_type dim0;
        size_type dim1;
        size_type dim2;
        size_type dim3;
        size_type dim4;
        size_type dim5;
        fftw_complex* data;
};



#endif // __ARRAY4D_H__
