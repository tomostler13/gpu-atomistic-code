// File: error.cpp
// Author:Tom Ostler
// Last-modified: 23 Oct 2015 17:45:51
#include "../inc/error.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
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
 /*   void errMessage(std::stringstream em)
    {
        std::string str=em.str();
        std::cerr << ". " << str.c_str() << "***\n" << std::endl;
        exit(0);
    }*/
    void errWarning(std::string em)
    {
        std::cerr << ". " << em << "***\n" << std::endl;
    }


}
