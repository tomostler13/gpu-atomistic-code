// File: error.h
// Author:Tom Ostler
// Last-modified: 23 Oct 2015 17:45:59
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#ifndef _ERROR_H_
#define _ERROR_H_
namespace error
{
    void errPreamble(std::string,int);
    void errMessage(std::string);
//    void errMessage(std::stringstream);
	void errWarning(std::string);

}
#endif /*_ERROR_H_*/
