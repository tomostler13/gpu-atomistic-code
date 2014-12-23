// File: exch.cpp
// Author: Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 15 Dec 2014 11:44:50
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/llg.h"
#include "../inc/matrix_conv.h"
#include "../inc/defines.h"
#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>
namespace exch
{
    unsigned int max_shells=0,diagnumdiag=0,offdiagnumdiag=0;
    Array<int> diagoffset,offdiagoffset;//for DIA neighbour list
    Array<int> checkdiag;
    Array<unsigned int> xadj,adjncy,offdiagxadj,offdiagadjncy;//for CSR neighbour list
    Array<double> dataxx,datayy,datazz,dataxz,dataxy,datayx,datayz,datazx,datazy;
    //evs = exchange vec scale
    Array<double> evs;
    Array3D<unsigned int> numint;
    Array2D<unsigned int> shell_list;
    Array4D<double> exchvec;
    Array5D<double> J;
    std::string enerType;
    bool outputJ,oem=false,rem=false;
    std::string readMethod,readFile,method;
    libconfig::Config exchcfg;
    void initExch(int argc,char *argv[])
    {

        readGlobalExch(argc,argv);
        if(geom::nspins>1)//if not then really we don't want any exchange anyway.
        {
            if(method=="permute" && rem==false)
            {
                permute(argc,argv);
            }//end of if(method=="permute") statement
            else if(method=="direct" && rem==false)
            {
                direct(argc,argv);
            }
            else if(rem)//then the interaction matrix has already been calculated and we just have to read it
            {
                if(config::exchm==2)
                {
                    std::ifstream ipem("csr_exch_mat.dat");
                    if(!ipem.is_open())
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("You requested to read in a CSR exchange matrix but it does not exist. Check it is called csr_exch_mag.dat");
                    }
                    unsigned int lns=0,adjsize;
                    ipem >> lns;
                    ipem >> adjsize;
                    if(lns!=geom::nspins)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("The number of spins (first line) in the CSR input file is not consistent with the system size you have generated. (First check failed)");
                    }
                    if(adjsize < 1)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("The size of adjncy (neighbour list) is zero. (Second check failed)");
                    }

                    FIXOUT(config::Info,"Size of xadj read as:" << lns+1 << std::endl);

                    FIXOUT(config::Info,"Resizing and filling xadj array:" << std::flush);
                    xadj.resize(geom::nspins+1);
                    xadj.IFill(0);
                    //read the xadj array
                    for(unsigned int i = 0 ; i < geom::nspins ; i++)
                    {
                        ipem >> xadj[i];
                    }
                    ipem >> xadj[geom::nspins];
                    SUCCESS(config::Info);
                    FIXOUT(config::Info,"Size of adjncy array read as:" << adjsize << std::endl);
                    FIXOUT(config::Info,"Resizing and filling adjncy array:" << std::flush);
                    adjncy.resize(adjsize);adjncy.IFill(0);
                    dataxx.resize(adjsize);dataxx.IFill(0);
                    datayy.resize(adjsize);datayy.IFill(0);
                    datazz.resize(adjsize);datazz.IFill(0);
                    for(unsigned int i = 0 ; i < geom::nspins ; i++)
                    {
                        for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
                        {
                            if(j > adjncy.size())
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("Stopping a seg fault. There is an error in the xadj lookup (adjncy would seg fault). (Third check failed)");
                            }
                            ipem >> adjncy[j] >> dataxx[j] >> datayy[j] >> datazz[j];
                        }
                    }
                    ipem.close();
                    if(ipem.is_open())
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errWarning("Could not close the csr_exch_mat.dat file.");
                    }
                }
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Exchange method not recognized");
            }
        }
    }
}
