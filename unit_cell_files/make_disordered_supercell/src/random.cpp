// File: random.cpp
// Author: Joe Barker
// Last-modified: 22 Nov 2014 21:28:30
#include <cmath>
#include "../../../inc/random.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
namespace Random{

  const double _rmax=1.0/4294967295.0;
  bool   _norm_sw=false;
  double _norm1,_norm2;
  unsigned int _x,_y;
  bool ri=false;

  inline unsigned int MWC1616()
  {
    assert(ri);
    _x=(18000*(_x&0xFFFF)+(_x>>16));
    _y=(30903*(_y&0xFFFF)+(_y>>16));
    return (_x<<16)+(_y&0xFFFF);
  }

  void seed(const unsigned int a,const unsigned int b)
  {
    if(a == 0 || b == 0){
		std::cerr << "Cannot use 0 seed for random number generator." << std::endl;
		exit(0);
    }
    _x = a;
    _y = b;
    ri=true;
  }


  double rand()
  {
    return MWC1616()*_rmax;
  }

  double normal()
  {
    double w,x1,x2;
    if(_norm_sw==false){
      for(;;){
        const double u=1.0-MWC1616()*_rmax;
        const double v=1.0-MWC1616()*_rmax;
          x1=2.0*u-1.0;
          x2=2.0*v-1.0;
          w=x1*x1+x2*x2;
          if(w<1.0) break;
      }

      w=sqrt((-2.0*log(w))/w);
      _norm1=w*x1;
      _norm2=w*x2;
      _norm_sw=true;
      return _norm1;
    } else{
      _norm_sw=false;
      return _norm2;
    }
  }
}
