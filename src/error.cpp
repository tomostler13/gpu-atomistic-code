// File: error.cpp
// Author:Tom Ostler
// Last-modified: 02 Dec 2012 20:06:08
#include "../inc/error.h"
#include <iostream>
#include <cstdlib>
#include <string>
namespace error
{
    void errPreamble(std::string f,int l)
    {
        std::cerr << "\n***Error in File: " << f << ", at line: " << l <<std::flush;
    }
    void errMessage(std::string em)
    {
        std::cerr << ". " << em << "***\n" << std::endl;
        exit(0);
    }
    void errWarning(std::string em)
    {
        std::cerr << ". " << em << "***\n" << std::endl;
    }


}
