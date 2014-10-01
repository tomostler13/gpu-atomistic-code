// File: array5d.h
// Author:Joe Barker
// Last-modified: 01 Oct 2014 18:39:26
#ifndef __ARRAY7D_H__
#define __ARRAY7D_H__

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
class Array7D
{
    public:
        typedef unsigned int size_type;
        typedef _Tp value_type;
        typedef _Tp* iterator;
        typedef const _Tp* const_iterator;
        typedef ptrdiff_t difference_type;

        Array7D() : dim0(0), dim1(0), dim2(0), dim3(0), dim4(0), dim5(0), dim6(0), data(0) {}

        Array7D(size_type d0, size_type d1, size_type d2, size_type d3, size_type d4, size_type d5, size_type d6)
            : dim0(d0), dim1(d1), dim2(d2), dim3(d3), dim4(d4), dim5(d5), dim6(d6), data(d0*d1*d2*d3*d4*d5*d6) {}

        ~Array7D(){data.clear();}

        inline void clear() {
            dim0 = 0;
            dim1 = 0;
            dim2 = 0;
            dim3 = 0;
            dim4 = 0;
            dim5 = 0;
            dim6 = 0;
            data.clear();
        }

        inline int getarrayelement(unsigned int i,unsigned int j,unsigned int k,unsigned int l,unsigned int m,unsigned int n,unsigned int o)
        {
            return((((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o);
        }
        inline int size()
        {
            return(dim0*dim1*dim2*dim3*dim4*dim5*dim6);
        }

        inline void resize(const size_type d0, const size_type d1, const size_type d2, const size_type d3, const size_type d4,const size_type d5, const size_type d6) {
            dim0 = d0; dim1 = d1; dim2 = d2; dim3 = d3; dim4 = d4; dim5=d5; dim6=d6;
            data.resize(d0*d1*d2*d3*d4*d5*d6);
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
                                    for(unsigned int o = 0 ; o < dim6 ; o++)
                                    {
                                        data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o]=a;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        inline _Tp& RESTRICT operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n, const size_type o) {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            assert( o < dim6 );
            return data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o];
        }

        inline const _Tp& operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n,const size_type o) const {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            assert( o < dim6 );
            return data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o];
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
        size_type dim6;
        std::vector<_Tp> data;
};

template <>
class Array7D<fftw_complex>
{
    public:
        typedef unsigned int size_type;
        typedef fftw_complex value_type;
        typedef fftw_complex* iterator;
        typedef const fftw_complex* const_iterator;
        typedef ptrdiff_t difference_type;

        Array7D() : dim0(0), dim1(0), dim2(0), dim3(0), dim4(0), dim5(0), dim6(0), data(0) {}

        Array7D(size_type d0, size_type d1, size_type d2, size_type d3, size_type d4, size_type d5, size_type d6)
        {resize(d0,d1,d2,d3,d4,d5,d6);}

        ~Array7D(){clear();}

        inline void clear() {
            dim0 = 0;
            dim1 = 0;
            dim2 = 0;
            dim3 = 0;
            dim4 = 0;
            dim5 = 0;
            dim6 = 0;
            fftw_free(data);
            data = NULL;
        }
        inline int size()
        {
            return(dim0*dim1*dim2*dim3*dim4*dim5*dim6);
        }
        inline void resize(const size_type d0, const size_type d1, const size_type d2, const size_type d3, const size_type d4,const size_type d5,const size_type d6) {
            dim0 = d0; dim1 = d1; dim2 = d2; dim3 = d3; dim4 = d4; dim5 = d5;dim6=d6;
            if(data != NULL){
                fftw_free(data);
            }
            data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*d0*d1*d2*d3*d4*d5*d6);
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
                                    for(unsigned int o = 0 ; o < dim6 ; o++)
                                    {
                                        data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o][0]=a;
                                        data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o][1]=a;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        inline fftw_complex& RESTRICT operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n,const size_type o) {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            assert( o < dim6 );
            return data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o];
        }

        inline const fftw_complex& operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n,const size_type o) const {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            assert( o < dim6 );
            return data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o];
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
        size_type dim6;
        fftw_complex* data;
};

template <>
class Array7D<fftwf_complex>
{
    public:
        typedef unsigned int size_type;
        typedef fftwf_complex value_type;
        typedef fftwf_complex* iterator;
        typedef const fftwf_complex* const_iterator;
        typedef ptrdiff_t difference_type;

        Array7D() : dim0(0), dim1(0), dim2(0), dim3(0), dim4(0), dim5(0), dim6(0), data(0) {}

        Array7D(size_type d0, size_type d1, size_type d2, size_type d3, size_type d4, size_type d5, size_type d6)
        {resize(d0,d1,d2,d3,d4,d5,d6);}

        ~Array7D(){clear();}

        inline void clear() {
            dim0 = 0;
            dim1 = 0;
            dim2 = 0;
            dim3 = 0;
            dim4 = 0;
            dim5 = 0;
            dim6 = 0;
            fftwf_free(data);
            data = NULL;
        }
        inline int size()
        {
            return(dim0*dim1*dim2*dim3*dim4*dim5*dim6);
        }
        inline void resize(const size_type d0, const size_type d1, const size_type d2, const size_type d3, const size_type d4,const size_type d5,const size_type d6) {
            dim0 = d0; dim1 = d1; dim2 = d2; dim3 = d3; dim4 = d4; dim5 = d5;dim6=d6;
            if(data != NULL){
                fftwf_free(data);
            }
            data = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*d0*d1*d2*d3*d4*d5*d6);
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
                                    for(unsigned int o = 0 ; o < dim6 ; o++)
                                    {
                                        data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o][0]=a;
                                        data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o][1]=a;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        inline fftwf_complex& RESTRICT operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n,const size_type o) {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            assert( o < dim6 );
            return data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o];
        }

        inline const fftwf_complex& operator()(const size_type i, const size_type j,
                const size_type k, const size_type l, const size_type m, const size_type n,const size_type o) const {
            assert( i < dim0 );
            assert( j < dim1 );
            assert( k < dim2 );
            assert( l < dim3 );
            assert( m < dim4 );
            assert( n < dim5 );
            assert( o < dim6 );
            return data[(((((i*dim1+j)*dim2+k)*dim3+l)*dim4+m)*dim5+n)*dim6+o];
        }

        inline fftwf_complex* __restrict__ ptr() {
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
        size_type dim6;
        fftwf_complex* data;
};


#endif // __ARRAY4D_H__
