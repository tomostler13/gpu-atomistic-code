// File: main.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 15 Jan 2013 19:34:25
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include "../inc/error.h"
#include "../inc/config.h"
int main(int argc,char *argv[])
{
    config::initConfig(argc,argv);
    return(0);
}
