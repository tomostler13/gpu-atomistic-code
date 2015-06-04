// File: exch_read_4_spin.cpp
// Author: Tom Ostler
// Created: 03 June 2015
// Last-modified: 04 Jun 2015 12:56:31
// If the 4 spin terms are included the routines
// in this source file will be called. It reads
// permutations of the quartets. It then mallocs
// the memory for containing the 4D sparse matrix
// interaction.
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
    void read4spin()
    {
        fsq.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_4s,3,3);
        fsqi.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_4s,3,3);
        numquart.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS());
        fsq.IFill(0);
        fsq.IFill(0);
        numquart.IFill(0);
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
            for(unsigned int j = 0 ; j < geom::ucm.GetNMS() ; j++)
            {
                std::stringstream sstr;
                sstr << "FourSpinExchange_" << i << "_" << j;
                std::string str=sstr.str();
                bool c4s=exchcfg.exists(str.c_str());

                if(!c4s)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting " << sstr << " does not exist. Check your config file.";
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
                libconfig::Setting &setting = exchcfg.lookup(str.c_str());
                config::printline(config::Info);
                if(setting.lookupValue("NumPermutations",numquart(i,j)))
                {
                    std::stringstream sstr;
                    sstr << "The number of four spin quartets between species " << i << " and " << j << " is:";
                    std::string str=sstr.str();
                    FIXOUT(config::Info,str << numquart(i,j) << std::endl);
                }
                //now read the quartets
                for(unsigned int q = 0 ; q < numquart(i,j) ; q++)
                {
                    for(unsigned int jkl = 0 ; jkl < 3 ; jkl++)
                    {
                        std::stringstream sstr;
                        sstr << "Permute" << q+1 << "_" << jkl+1;
                        std::string str=sstr.str();
                        for(unsigned int c = 0 ; c < 3 ; c++)
                        {
                            try
                            {
                                fsq(i,j,q,jkl,c)=setting[str.c_str()][c];
                                fsqi(i,j,q,jkl,c)=static_cast<int>(fsq(i,j,q,jkl,c)*static_cast<double>(geom::Nk[c]*evs[c]+0.5));
                            }
                            catch(const libconfig::SettingNotFoundException &snf)
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                std::stringstream nss;
                                nss << "Setting not found exception caught. Setting " << snf.getPath();
                                std::string ns=nss.str();
                                error::errMessage(ns);
                            }
                        }
                        std::stringstream qss;
                        qss << "Quartet " << q << " (number in quartet out of 3) " << jkl+1;
                        std::string qs=qss.str();
                        FIXOUTVEC(config::Info,qs,fsq(i,j,q,jkl,0),fsq(i,j,q,jkl,1),fsq(i,j,q,jkl,2));
                        FIXOUTVEC(config::Info,"On integer mesh:",fsqi(i,j,q,jkl,0),fsqi(i,j,q,jkl,1),fsqi(i,j,q,jkl,2));


                    }
                    config::printline(config::Info);
                }

            }
        }
    }

    void setup4SpinCSR()
    {
        xadj_j.resize(geom::nspins+1);
        xadj_k.resize(geom::nspins+1);
        xadj_l.resize(geom::nspins+1);
        std::vector<unsigned int> tempadjncy_j,tempadjncy_k,tempadjncy_l;
        //loop over all spins and find the neighbours
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {

        }
    }
}
