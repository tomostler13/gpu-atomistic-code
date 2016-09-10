// File: error.cpp
// Author:Tom Ostler
// Last-modified: 25 Feb 2016 21:33:28
#include "../../inc/error.h"
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
