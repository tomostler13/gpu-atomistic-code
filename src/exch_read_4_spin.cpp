// File: exch_read_4_spin.cpp
// Author: Tom Ostler
// Created: 03 June 2015
// Last-modified: 05 Jun 2015 20:01:55
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
#include "../inc/fields.h"
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
                                if(fsq(i,j,q,jkl,c)<0)
                                {
                                    fsqi(i,j,q,jkl,c)=static_cast<int>(fsq(i,j,q,jkl,c)*static_cast<double>(geom::Nk[c]*evs[c]-0.5));
                                }
                                else
                                {
                                     fsqi(i,j,q,jkl,c)=static_cast<int>(fsq(i,j,q,jkl,c)*static_cast<double>(geom::Nk[c]*evs[c]+0.5));
                                }
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
        setup4SpinCSR();
    }

    void setup4SpinCSR()
    {
        xadj_j.resize(geom::nspins+1);
        xadj_k.resize(geom::nspins+1);
        xadj_l.resize(geom::nspins+1);
        std::vector<unsigned int> tempadjncy_j,tempadjncy_k,tempadjncy_l;
        //loop over all spins and find the neighbours
        unsigned int xadj_counter=0,adjncy_count_j=0,adjncy_count_k=0,adjncy_count_l=0;
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            xadj_j[xadj_counter]=adjncy_count_j;
            xadj_k[xadj_counter]=adjncy_count_k;
            xadj_l[xadj_counter]=adjncy_count_l;
            unsigned int sl=geom::sublattice(i);
            //get the atom in the unit cell
            unsigned int aiuc = geom::lu(i,4);

            //store i's position
            int mypos[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
            for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
            {
                for(unsigned int q = 0 ; q < numquart(sl,s1) ; q++)
                {
                    //this could be tidied up by doing a loop
                    //over j, k and l
                    int lookupj[3]={fsqi(sl,s1,q,0,0),fsqi(sl,s1,q,0,1),fsqi(sl,s1,q,0,2)};
                    int lookupk[3]={fsqi(sl,s1,q,1,0),fsqi(sl,s1,q,1,1),fsqi(sl,s1,q,1,2)};
                    int lookupl[3]={fsqi(sl,s1,q,2,0),fsqi(sl,s1,q,2,1),fsqi(sl,s1,q,2,2)};
                    int lookupvecj[3]={mypos[0]+lookupj[0],mypos[1]+lookupj[1],mypos[2]+lookupj[2]};
                    int lookupveck[3]={mypos[0]+lookupk[0],mypos[1]+lookupk[1],mypos[2]+lookupk[2]};
                    int lookupvecl[3]={mypos[0]+lookupl[0],mypos[1]+lookupl[1],mypos[2]+lookupl[2]};


                    unsigned int check_lookupj=0,check_lookupk=0,check_lookupl=0;
                    //check pbc's for j
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        if(lookupvecj[xyz]<0 && config::pbc[xyz]==true)
                        {
                            lookupvecj[xyz]=geom::dim[xyz]*geom::Nk[xyz]+lookupvecj[xyz];
                            check_lookupj++;
                        }
                        else if(lookupvecj[xyz]>=geom::dim[xyz]*geom::Nk[xyz] && config::pbc[xyz]==true)
                        {
                            lookupvecj[xyz]=lookupvecj[xyz]-geom::dim[xyz]*geom::Nk[xyz];
                            check_lookupj++;
                        }
                        else if(lookupvecj[xyz]>=0 && lookupvecj[xyz]<geom::dim[xyz]*geom::Nk[xyz])
                        {
                            check_lookupj++;
                        }
                    }//end of xyz loop for j
                    //check pbc's for k
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        if(lookupveck[xyz]<0 && config::pbc[xyz]==true)
                        {
                            lookupveck[xyz]=geom::dim[xyz]*geom::Nk[xyz]+lookupveck[xyz];
                            check_lookupk++;
                        }
                        else if(lookupveck[xyz]>=geom::dim[xyz]*geom::Nk[xyz] && config::pbc[xyz]==true)
                        {
                            lookupveck[xyz]=lookupveck[xyz]-geom::dim[xyz]*geom::Nk[xyz];
                            check_lookupk++;
                        }
                        else if(lookupveck[xyz]>=0 && lookupveck[xyz]<geom::dim[xyz]*geom::Nk[xyz])
                        {
                            check_lookupk++;
                        }
                    }//end of xyz loop for k
                    //check pbc's for l
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        if(lookupvecl[xyz]<0 && config::pbc[xyz]==true)
                        {
                            lookupvecl[xyz]=geom::dim[xyz]*geom::Nk[xyz]+lookupvecl[xyz];
                            check_lookupl++;
                        }
                        else if(lookupvecl[xyz]>=geom::dim[xyz]*geom::Nk[xyz] && config::pbc[xyz]==true)
                        {
                            lookupvecl[xyz]=lookupvecl[xyz]-geom::dim[xyz]*geom::Nk[xyz];
                            check_lookupl++;
                        }
                        else if(lookupvecl[xyz]>=0 && lookupvecl[xyz]<geom::dim[xyz]*geom::Nk[xyz])
                        {
                            check_lookupl++;
                        }
                    }//end of xyz loop for l
                    if(check_lookupj==3)
                    {
                        unsigned int spec=geom::coords(lookupvecj[0],lookupvecj[1],lookupvecj[2],1);
                        if(spec==sl)
                        {
                            unsigned int neighj=geom::coords(lookupvecj[0],lookupvecj[1],lookupvecj[2],0);
                            tempadjncy_j.push_back(neighj);
                            tempadjncy_j[adjncy_count_j]=neighj;
                            adjncy_count_j++;
                        }//end of check that atom exists and is correct spec
                    }
                    if(check_lookupk==3)
                    {
                        unsigned int spec=geom::coords(lookupveck[0],lookupveck[1],lookupveck[2],1);
                        if(spec==sl)
                        {
                            unsigned int neighk=geom::coords(lookupveck[0],lookupveck[1],lookupveck[2],0);
                            tempadjncy_k.push_back(neighk);
                            tempadjncy_k[adjncy_count_k]=neighk;
                            adjncy_count_k++;
                        }//end of check that atom exists and is correct spec
                    }
                    if(check_lookupl==3)
                    {
                        unsigned int spec=geom::coords(lookupvecl[0],lookupvecl[1],lookupvecl[2],1);
                        if(spec==sl)
                        {
                            unsigned int neighl=geom::coords(lookupvecl[0],lookupvecl[1],lookupvecl[2],0);
                            tempadjncy_l.push_back(neighl);
                            tempadjncy_l[adjncy_count_l]=neighl;
                            adjncy_count_l++;
                        }//end of check that atom exists and is correct spec
                    }
                }//end of q loop
            }//end of s1 loop
            xadj_counter++;


        }//end of i loop
        xadj_j[geom::nspins]=adjncy_count_j;
        xadj_k[geom::nspins]=adjncy_count_k;
        xadj_l[geom::nspins]=adjncy_count_l;

        adjncy_j.resize(adjncy_count_j);
        adjncy_k.resize(adjncy_count_k);
        adjncy_l.resize(adjncy_count_l);
        for(unsigned int i = 0 ; i < adjncy_count_j ; i++)
        {
            adjncy_j[i]=tempadjncy_j[i];
        }
        for(unsigned int i = 0 ; i < adjncy_count_k ; i++)
        {
            adjncy_k[i]=tempadjncy_k[i];
        }
        for(unsigned int i = 0 ; i < adjncy_count_l ; i++)
        {
            adjncy_l[i]=tempadjncy_l[i];
        }

        FIXOUT(config::Info,"Four Spin -> size of adjncy_j:" << adjncy_j.size() << std::endl);
        FIXOUT(config::Info,"Four Spin -> size of adjncy_k:" << adjncy_k.size() << std::endl);
        FIXOUT(config::Info,"Four Spin -> size of adjncy_l:" << adjncy_l.size() << std::endl);
    }//end of void setup4spinCSR function
}
