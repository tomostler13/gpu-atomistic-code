// File: error.h
// Author:Tom Ostler
// Last-modified: 02 Dec 2012 20:08:36
#include <string>
#include <iostream>
#include <cstdlib>
#ifndef _ERROR_H_
#define _ERROR_H_
namespace error
{
    void errPreamble(std::string,int);
    void errMessage(std::string);
	void errWarning(std::string);

}
#endif /*_ERROR_H_*/
