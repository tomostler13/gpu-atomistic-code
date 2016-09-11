// File: array4d.h
// Author:Joe Barker
// Last-modified: 23 Nov 2012 17:10:05

#ifndef __ARRAY4D_H__
#define __ARRAY4D_H__

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
class Array4D
{
  public:
    typedef unsigned int size_type;
    typedef _Tp value_type;
    typedef _Tp* iterator;
    typedef const _Tp* const_iterator;
    typedef ptrdiff_t difference_type;

    Array4D() : dim0(0), dim1(0), dim2(0), dim3(0), data(0) {}

    Array4D(size_type d0, size_type d1, size_type d2, size_type d3)
      : dim0(d0), dim1(d1), dim2(d2), dim3(d3), data(d0*d1*d2*d3) {}

    ~Array4D(){data.clear();}

    inline void clear() {
      dim0 = 0;
      dim1 = 0;
      dim2 = 0;
      dim3 = 0;
      data.clear();
    }

    inline void resize(const size_type d0, const size_type d1, const size_type d2, const size_type d3) {
      dim0 = d0; dim1 = d1; dim2 = d2; dim3 = d3;
      data.resize(d0*d1*d2*d3);
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
                        data[((i*dim1+j)*dim2+k)*dim3+l]=a;
                    }
                }
            }
        }
    }


    inline _Tp& RESTRICT operator()(const size_type i, const size_type j,
        const size_type k, const size_type l) {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      assert( l < dim3 );
      return data[((i*dim1+j)*dim2+k)*dim3+l];
    }
    inline _Tp* RESTRICT ptr() {
      return &(data[0]);
    }

    inline const _Tp& operator()(const size_type i, const size_type j,
        const size_type k, const size_type l) const {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      assert( l < dim3 );
      return data[((i*dim1+j)*dim2+k)*dim3+l];
    }

//    inline size_type size(const size_type i) const { return dim[i]; }

  private:
    size_type dim0;
    size_type dim1;
    size_type dim2;
    size_type dim3;
    std::vector<_Tp> data;
};

template <>
class Array4D<fftw_complex>
{
  public:
    typedef unsigned int size_type;
    typedef fftw_complex value_type;
    typedef fftw_complex* iterator;
    typedef const fftw_complex* const_iterator;
    typedef ptrdiff_t difference_type;

    Array4D() : dim0(0), dim1(0), dim2(0), dim3(0), data(0) {}

    Array4D(size_type d0, size_type d1, size_type d2, size_type d3)
      {resize(d0,d1,d2,d3);}

    ~Array4D(){clear();}

    inline void clear() {
      dim0 = 0;
      dim1 = 0;
      dim2 = 0;
      dim3 = 0;
      fftw_free(data);
      data = NULL;
    }

    inline void resize(const size_type d0, const size_type d1, const size_type d2, const size_type d3) {
      dim0 = d0; dim1 = d1; dim2 = d2; dim3 = d3;
      if(data != NULL){
        fftw_free(data);
      }
      data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*d0*d1*d2*d3);
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
                        data[((i*dim1+j)*dim2+k)*dim3+l][0]=a;
                        data[((i*dim1+j)*dim2+k)*dim3+l][1]=a;
                    }
                }
            }
        }
    }

    inline fftw_complex& RESTRICT operator()(const size_type i, const size_type j,
        const size_type k, const size_type l) {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      assert( l < dim3 );
      return data[((i*dim1+j)*dim2+k)*dim3+l];
    }

    inline const fftw_complex& operator()(const size_type i, const size_type j,
        const size_type k, const size_type l) const {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      assert( l < dim3 );
      return data[((i*dim1+j)*dim2+k)*dim3+l];
    }

	inline fftw_complex* __restrict__ ptr()
	{
		return data;
	}
//    inline size_type size(const size_type i) const { return dim[i]; }

  private:
    size_type dim0;
    size_type dim1;
    size_type dim2;
    size_type dim3;
    fftw_complex* data;
};



#endif // __ARRAY4D_H__
