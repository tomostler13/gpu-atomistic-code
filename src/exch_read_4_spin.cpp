// File: exch_read_4_spin.cpp
// Author: Tom Ostler
// Created: 03 June 2015
// Last-modified: 08 Dec 2015 15:53:31
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
        fsq.resize(geom::ucm.GetNMS(),max_4s,3,3);
        fsqi.resize(geom::ucm.GetNMS(),max_4s,3,3);
        numquart.resize(geom::ucm.GetNMS());
        JQ.resize(geom::ucm.GetNMS());
        //At the moment there can only be one interface per species
        Interface.resize(geom::ucm.GetNMS());
        fsq.IFill(0);
        fsq.IFill(0);
        numquart.IFill(0);
        JQ.IFill(0);
        Interface.IFill(-1);
        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
        {
                std::stringstream sstr;
                sstr << "FourSpinExchange_" << i;
                std::string str=sstr.str();
                bool c4s=exchcfg.exists(str.c_str());

                if(!c4s)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting " << str << " does not exist. Check your config file.";
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
                libconfig::Setting &setting = exchcfg.lookup(str.c_str());
                config::printline(config::Info);
                if(setting.lookupValue("NumPermutations",numquart(i)))
                {
                    std::stringstream sstr;
                    sstr << "The number of four spin quartets between species " << i << " is:";
                    std::string str=sstr.str();
                    FIXOUT(config::Info,str << numquart(i) << std::endl);
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream sstr;
                    sstr << "Could not read the number of quartets between species " << i;
                    std::string str=sstr.str();
                    error::errMessage(str);
                }
                if(setting.lookupValue("JQ",JQ(i)))
                {
                    std::stringstream sstr;
                    sstr << "The four spin exchange between species " << i << " is:";
                    std::string str=sstr.str();
                    FIXOUT(config::Info,str << JQ(i) << std::endl);
                    JQ(i)=JQ(i)/(geom::ucm.GetMu(i)*9.27e-24);
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream sstr;
                    sstr << "Could not read the four spin exchange between species " << i;
                    std::string str=sstr.str();
                    error::errMessage(str);
                }
                if(setting.lookupValue("InterfaceSpecies",Interface[i]))
                {
                    FIXOUT(config::Info,"Interface with species:" << Interface[i] << std::endl);
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not read sublattice for interfaces");
                }
                //now read the quartets
                for(unsigned int q = 0 ; q < numquart(i) ; q++)
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
                                fsq(i,q,jkl,c)=setting[str.c_str()][c];
/*                                if(fsq(i,j,q,jkl,c)<0)
                                {
                                    fsqi(i,j,q,jkl,c)=static_cast<int>(fsq(i,j,q,jkl,c)*static_cast<double>(geom::Nk[c]*evs[c]-0.5));
                                }
                                else
                                {*/
                                     fsqi(i,q,jkl,c)=static_cast<int>(fsq(i,q,jkl,c)*static_cast<double>(geom::Nk[c]*evs[c]+0.5));
                                //}
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
                        FIXOUTVEC(config::Info,qs,fsq(i,q,jkl,0),fsq(i,q,jkl,1),fsq(i,q,jkl,2));
                        FIXOUTVEC(config::Info,"On integer mesh:",fsqi(i,q,jkl,0),fsqi(i,q,jkl,1),fsqi(i,q,jkl,2));


                    }
                    config::printline(config::Info);

            }
        }
        setup4SpinCSR();
    }

    void setup4SpinCSR()
    {
        xadj_j.resize(geom::nspins+1);
        std::vector<unsigned int> tempadjncy_j,tempadjncy_k,tempadjncy_l;
        //loop over all spins and find the neighbours
        unsigned int xadj_counter=0,adjncy_count_j=0,adjncy_count_k=0,adjncy_count_l=0;
        unsigned int totquartets=0,missingquartets=0;

        Array2D<unsigned int> checkq;
        checkq.resize(geom::nspins,max_4s);
        checkq.IFill(0);
        Array<unsigned int> countq;
        countq.resize(geom::nspins);
        countq.IFill(0);
        if(!rem)
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                xadj_j[xadj_counter]=adjncy_count_j;
                unsigned int sl=geom::sublattice(i);

                //store i's position
                int mypos[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};

                for(unsigned int q = 0 ; q < numquart(sl) ; q++)
                {
                    //this could be tidied up by doing a loop
                    //over j, k and l
                    int lookupj[3]={fsqi(sl,q,0,0),fsqi(sl,q,0,1),fsqi(sl,q,0,2)};
                    int lookupk[3]={fsqi(sl,q,1,0),fsqi(sl,q,1,1),fsqi(sl,q,1,2)};
                    int lookupl[3]={fsqi(sl,q,2,0),fsqi(sl,q,2,1),fsqi(sl,q,2,2)};
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
                    if((check_lookupj+check_lookupk+check_lookupl)==9)//if there aren't four spins in the quartet we ignore the entire 4 spin term
                    {
                        unsigned int specj=geom::coords(lookupvecj[0],lookupvecj[1],lookupvecj[2],1);
                        unsigned int speck=geom::coords(lookupveck[0],lookupveck[1],lookupveck[2],1);
                        unsigned int specl=geom::coords(lookupvecl[0],lookupvecl[1],lookupvecl[2],1);
                        if((specj==sl || specj==Interface[sl]) && (speck==sl || speck==Interface[sl]) && (specl==sl || specl==Interface[sl]) && checkq(i,q)==0)
                        {
                            unsigned int neighj=geom::coords(lookupvecj[0],lookupvecj[1],lookupvecj[2],0);
                            tempadjncy_j.push_back(neighj);
                            tempadjncy_j[adjncy_count_j]=neighj;
                            adjncy_count_j++;
                            unsigned int neighk=geom::coords(lookupveck[0],lookupveck[1],lookupveck[2],0);
                            tempadjncy_k.push_back(neighk);
                            tempadjncy_k[adjncy_count_k]=neighk;
                            adjncy_count_k++;
                            unsigned int neighl=geom::coords(lookupvecl[0],lookupvecl[1],lookupvecl[2],0);
                            tempadjncy_l.push_back(neighl);
                            tempadjncy_l[adjncy_count_l]=neighl;
                            adjncy_count_l++;
                            checkq(i,q)=1;
                            countq[i]++;

                            totquartets++;
                        }
                        else
                        {
                            missingquartets++;
                        }
                    }//end of check of sum of check_lookup equal to 9
                }//end of q loop
                xadj_counter++;

                config::Log << "spin " << i << " has " << countq[i] << " quartets" << std::endl;
            }//end of i loop
            xadj_j[geom::nspins]=adjncy_count_j;
            FIXOUT(config::Info,"Total number of quartets:" << totquartets << std::endl);
            FIXOUT(config::Info,"Total number of missing quartets:" << missingquartets << std::endl);
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
        }
        else
        {
            std::ifstream ipem("4spin_csr_exch_mat.dat");
            if(!ipem.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You requested to read in a CSR exchange matrix but it does not exist. Check it is called 4spin_csr_exch_mat.dat");
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

            FIXOUT(config::Info,"4spin -> Size of xadj read as:" << lns+1 << std::endl);

            FIXOUT(config::Info,"4spin -> Resizing and filling xadj array:" << std::flush);
            xadj_j.resize(geom::nspins+1);
            xadj_j.IFill(0);
            //read the xadj array
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                ipem >> xadj_j[i];
            }
            ipem >> xadj_j[geom::nspins];
            SUCCESS(config::Info);
            FIXOUT(config::Info,"4spin -> Size of adjncy array read as:" << adjsize << std::endl);
            FIXOUT(config::Info,"4spin -> Resizing and filling adjncy arrays:" << std::flush);
            adjncy_j.resize(adjsize);adjncy_j.IFill(0);
            adjncy_k.resize(adjsize);adjncy_k.IFill(0);
            adjncy_l.resize(adjsize);adjncy_l.IFill(0);
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                for(unsigned int j = xadj_j[i] ; j < xadj_j[i+1] ; j++)
                {
                    if(j > adjncy_j.size())
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("Stopping a seg fault. There is an error in the xadj lookup (adjncy would seg fault). (Third check failed)");
                    }
                    ipem >> adjncy_j[j] >> adjncy_k[j] >> adjncy_l[j];
                }
            }
            ipem.close();
            if(ipem.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errWarning("Could not close the 4spin_csr_exch_mat.dat file.");
            }
        }
        FIXOUT(config::Info,"Four Spin -> size of adjncy_j:" << adjncy_j.size() << std::endl);
        FIXOUT(config::Info,"Four Spin -> size of adjncy_k:" << adjncy_k.size() << std::endl);
        FIXOUT(config::Info,"Four Spin -> size of adjncy_l:" << adjncy_l.size() << std::endl);
        if(oem)
        {
            FIXOUT(config::Info,"Outputting 4 spin exchange matrix:" << std::flush);
            std::ofstream opem("4spin_csr_exch_mat.dat");
            if(!opem.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open file for outputting the 4spin csr matrix.");
            }
            opem << geom::nspins << std::endl;
            opem << adjncy_j.size() << std::endl;
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                opem << xadj_j[i] << "\t";
            }
            opem << xadj_j[geom::nspins] << std::endl;
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                for(unsigned int j = xadj_j[i] ; j < xadj_j[i+1] ; j++)
                {
                    opem << adjncy_j[j] << "\t" << adjncy_k[j] << "\t" << adjncy_l[j] << "\t";
                }
            }
            opem.close();
            if(opem.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errWarning("Could not close the file for outputting the 4spin CSR exchange matrix.");
            }
            config::Info << "Done" << std::endl;
        }

    }//end of void setup4spinCSR function
}
