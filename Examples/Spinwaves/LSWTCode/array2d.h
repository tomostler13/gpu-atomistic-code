// File: array2d.h
// Author:Joe Barker
// Last-modified: 31 Dec 2012 12:10:25

#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

#include <vector>
#include <cassert>
#include <cstddef>
#include <fftw3.h>

// C ordered array
template <typename _Tp>
class Array2D
{
  public:
    typedef unsigned int size_type;
    typedef _Tp value_type;
    typedef _Tp* iterator;
    typedef const _Tp* const_iterator;
    typedef ptrdiff_t difference_type;

    Array2D() : dim0(0), dim1(0), data(0) {}

    Array2D(size_type d0, size_type d1)
      : dim0(d0), dim1(d1), data(d0*d1) {}

    ~Array2D(){data.clear();}

    inline void clear() {
      dim0 = 0;
      dim1 = 0;
      data.clear();
    }
   ///IFill means initialize or fill
    inline void IFill(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                data[i*dim1+j]=a;
            }
        }
    }

    inline void resize(size_type d0, size_type d1) {
      dim0 = d0; dim1 = d1;
      data.resize(d0*d1);
    }

    inline _Tp& RESTRICT operator()(const size_type i, const size_type j) {
      assert( i < dim0 );
      assert( j < dim1 );
      return data[i*dim1+j];
    }

    inline const _Tp& operator()(const size_type i, const size_type j) const {
      assert( i < dim0 );
      assert( j < dim1 );
      return data[i*dim1+j];
    }

    inline _Tp* RESTRICT ptr() {
      return &(data[0]);
    }

    inline size_type size() const { return (dim0*dim1); }

    private:
    size_type dim0;
    size_type dim1;
    std::vector<_Tp> data;
    };

    template <>
    class Array2D<fftw_complex>
    {
    public:
    typedef unsigned int size_type;
    typedef fftw_complex value_type;
    typedef fftw_complex* iterator;
    typedef const fftw_complex* const_iterator;
    typedef ptrdiff_t difference_type;

    Array2D() : dim0(0), dim1(0), data(NULL) {}

    Array2D(size_type d0, size_type d1)
      {resize(d0,d1);}

    ~Array2D(){clear();}

    inline void clear() {
      dim0 = 0;
      dim1 = 0;
      fftw_free(data);
      data = NULL;
    }
   ///IFill means initialize or fill
    inline void IFill(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                data[i*dim1+j][0]=a;
                data[i*dim1+j][1]=a;
            }
        }
    }

    inline void resize(size_type d0, size_type d1) {
      dim0 = d0; dim1 = d1;
      if(data != NULL){
          fftw_free(data);
      }
      data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*d0*d1);
    }

    inline fftw_complex& RESTRICT operator()(const size_type i, const size_type j) {
      assert( i < dim0 );
      assert( j < dim1 );
      return data[i*dim1+j];
    }

    inline const fftw_complex& operator()(const size_type i, const size_type j) const {
      assert( i < dim0 );
      assert( j < dim1 );
      return data[i*dim1+j];
    }

    inline fftw_complex* RESTRICT ptr() {
      return data;
    }

    //    inline size_type size(const size_type i) const { return dim[i]; }

    private:
    size_type dim0;
    size_type dim1;
    fftw_complex* data;
    };

    #endif // __ARRAY2D_H__
