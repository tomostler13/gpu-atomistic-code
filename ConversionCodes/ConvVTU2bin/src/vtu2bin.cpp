// File: vtu2bin.cpp
// Author:Tom Ostler
// Created: 22 June 2018
// Last-modified: 22 Jun 2018 08:58:11

//The purpose of this section of code is to take the output of the vtu
//files and write one file for the positions and then several files
//containing the relative spin components as binary
#include <cmath>
#include <iostream>
#include <libconfig.h++>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <fftw3.h>
#include "../../../inc/error.h"
#include "../../../inc/array.h"
#include "../../../inc/array2d.h"
#include "../../../inc/array3d.h"
#include "../../../inc/array4d.h"
#include "../../../inc/defines.h"
#include "../../../inc/random.h"
#include "../../inc/inputs.h"
int main(int argc,char *argv[])
{
    Array<double> sx,sy,sz,posx,posy,posz;
    std::cout << "#Reading config file..." << std::flush;
    inputs::readcff(argc,argv);

    std::cout << "Done" << std::endl;

    std::cout << "#Resizing magnetization and spin arrays..." << std::flush;
    sx.resize(inputs::nspins);
    sy.resize(inputs::nspins);
    sz.resize(inputs::nspins);

    sx.IFill(0);
    sy.IFill(0);
    sz.IFill(0);

    posx.resize(inputs::nspins);
    posy.resize(inputs::nspins);
    posz.resize(inputs::nspins);

    posx.IFill(0);
    posy.IFill(1);
    posz.IFill(2);


    std::cout << "Done" << std::endl;

    for(unsigned int t = inputs::lt ; t < inputs::ut+1 ; t+=inputs::dt)
    {
        std::stringstream sstr;
        sstr << inputs::fb << t << ".vtu";
        std::string str=sstr.str();

        std::cout << "#Opening file " << str << "..." << std::flush;
        std::ifstream ifs(str.c_str());
        if(!ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open input file");
        }
        std::cout << "done" << std::endl;
        std::string dumpline;
        //get rid of the top six lines
        for(unsigned int i = 0 ; i < 6 ; i++)
        {
            std::getline(ifs,dumpline);
//            std::cout << dumpline << std::endl;
        }
        //get the spin vectors (sx,sy,sz)
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            ifs >> sx[i] >> sy[i] >> sz[i];
        }
//        std::cout << spins(inputs::nspins-1,0) << std::endl;
        //get rid of the middle six lines
        for(unsigned int i = 0 ; i < 7 ; i++)
        {
            std::getline(ifs,dumpline);
//            std::cout << dumpline << std::endl;
        }
//        std::cin.get();
        //get the coordinates if we are looking at the first file
        if(t == inputs::lt)
        {
            for(unsigned int i = 0 ; i < inputs::nspins ; i++)
            {
                ifs >> posx[i] >> posy[i] >> posz[i];
            }
        }
        ifs.close();
        if(ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not close the ifs file");
        }
        std::stringstream osstr;
        osstr << "spins_" << t << ".bin";
        std::string ostr=osstr.str();
        std::ofstream ofs(ostr.c_str(),std::ios::out | std::ios::binary);
        if(!ofs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Couldn't open an output file");
        }
        else
        {
            ofs.write((char *) sx.ptr(),inputs::nspins*sizeof(double));
            ofs.write((char *) sy.ptr(),inputs::nspins*sizeof(double));
            ofs.write((char *) sz.ptr(),inputs::nspins*sizeof(double));
        }

        ofs.close();
    }
    std::ofstream posofs("pos.bin",std::ios::out | std::ios::binary);
    if(!posofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not open file for writing position, trying to write to cout");
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            std::cout << posx[i] << "\t" << posy[i] << "\t" << posz[i] << std::endl;
        }
    }
    else
    {
        posofs.write((char *) posx.ptr(),inputs::nspins*sizeof(double));
        posofs.write((char *) posy.ptr(),inputs::nspins*sizeof(double));
        posofs.write((char *) posz.ptr(),inputs::nspins*sizeof(double));
    }
    posofs.close();


    return(EXIT_SUCCESS);
}
