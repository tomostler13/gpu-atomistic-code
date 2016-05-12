// File: exch.cpp
// Author: Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 12 May 2016 18:32:19
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
    unsigned int max_shells=0,diagnumdiag=0,offdiagnumdiag=0,max_int=0,max_4s=0;
    Array<int> diagoffset,offdiagoffset;//for DIA neighbour list
    Array<int> checkdiag;
    Array<unsigned int> xadj,adjncy,offdiagxadj,offdiagadjncy,xadj_j,xadj_k,xadj_l,adjncy_j,adjncy_k,adjncy_l;//for CSR neighbour list
    Array<double> dataxx,datayy,datazz,dataxz,dataxy,datayx,datayz,datazx,datazy;
    //evs = exchange vec scale
    Array<double> evs;
    Array3D<unsigned int> numint;
    Array2D<unsigned int> shell_list;
    Array<unsigned int> numquart;
    Array<double> JQ;
    Array4D<double> exchvec;
    Array4D<int> fsqi;
    Array5D<double> J;
    Array4D<double> fsq;
    Array<int> Interface;
    std::string enerType,exchMatFN;
    bool outputJ,oem=false,rem=false,cutexch=false,inc4spin=false,eaem=false;
    //cut off of exchange in m
    double rcut=1.0;
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
            else if(method=="mapint" && rem==false)
            {
                mapint(argc,argv);
            }
            else if(rem)//then the interaction matrix has already been calculated and we just have to read it
            {
                if(config::exchm==2)
                {
                    std::ifstream ipem(exchMatFN.c_str());
                    if(!ipem.is_open())
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("You requested to read in a CSR exchange matrix but it does not exist. Check it is called csr_exch_mat.dat");
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

//                    std::ofstream cstr("check.mat");
//                    cstr << lns << std::endl;
//                    cstr << adjsize << std::endl;
                    //read the xadj array
                    for(unsigned int i = 0 ; i < geom::nspins ; i++)
                    {
                        ipem >> xadj[i];
//                        cstr << xadj[i] << "\t";
                    }
                    ipem >> xadj[geom::nspins];
//                    cstr << xadj[geom::nspins] << std::endl;
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
                            ipem >> std::setprecision(16) >> adjncy[j] >> dataxx[j] >> datayy[j] >> datazz[j];
                            //cstr << adjncy[j] << "\t" << dataxx[j] << "\t" << datayy[j] << "\t" << datazz[j] << "\t";
                        }
                    }
//                    cstr.close();
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
