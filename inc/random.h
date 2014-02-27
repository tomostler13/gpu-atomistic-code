// File: random.cpp
// Author: Joe Barker
// Last-modified: 05 Dec 2012 16:09:31
#ifndef _RANDOM_H_
#define _RANDOM_H_
namespace Random{
  extern bool ri;
  void seed(const unsigned int, const unsigned int);
  double rand();
  double normal();
}
#endif /*_RANDOM_H_*/
