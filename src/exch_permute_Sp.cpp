// File: exch_permute_Sp.cpp
// Author: Tom Ostler
// Created: 05 Dec 2014
// Last-modified: 15 Dec 2015 18:57:58
// This source file was added to tidy up the file exch.cpp
// because it was becoming cumbersome to work with. The
// routines here take the exchange shell list and permute
// the vectors and populate a sparse matrix.
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
    void exchpermute()
    {
        std::ofstream opJ;
        if(outputJ)
        {
            opJ.open("J.dat");
            if(!opJ.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open file J.dat for writing the 2D interaction matrix to.");
            }
        }
        //temporary arrays in format to store the adjncy as we don't know
        //how big it is before hand
        std::vector<unsigned int> tadjncy;
        std::vector< std::vector< std::vector< double > > > tdata;
        tdata.resize(3);
        tdata[0].resize(3);
        tdata[1].resize(3);
        tdata[2].resize(3);
        unsigned int adjcount=0;
        xadj.resize(geom::nspins+1);
        xadj.IFill(0);
        xadj[0]=0;
        unsigned int xadj_counter=0,adjncy_counter=0;

        //loop over the spins and find the neighbours
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            xadj[xadj_counter]=adjncy_counter;
            //what is my species (or sublattice)
            unsigned sl=geom::sublattice(i);
            //find my atom number in the unit cell so that
            //we can lookup the moment from the unit cell info
            unsigned int aiuc=geom::lu(i,4);
            //store atom i's position
            int mypos[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};

            Array3D<unsigned int> check;
            check.resize(geom::dim[0]*geom::Nk[0],geom::dim[1]*geom::Nk[1],geom::dim[2]*geom::Nk[2]);
            check.IFill(0);

            //this is a loop over the species.
            for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
            {
                for(unsigned int shell = 0 ; shell < shell_list(sl,s1) ; shell++)
                {
                    unsigned int lookup[3]={static_cast<unsigned int>(exchvec(sl,s1,shell,0)*evs[0]*static_cast<double>(geom::Nk[0])+0.5),
                        static_cast<unsigned int>(exchvec(sl,s1,shell,1)*evs[1]*static_cast<double>(geom::Nk[1])+0.5),
                        static_cast<unsigned int>(exchvec(sl,s1,shell,2)*evs[2]*static_cast<double>(geom::Nk[2])+0.5)};
//                        std::cout << "atom = " << i << " spec = " << s1 << " shell = " << shell << " lookup = " << lookup[0] << "\t" << lookup[1] << "\t" << lookup[2] << std::endl;
                    for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                    {
                        //reference array
                        int rc[3]={lookup[wrap%3],lookup[(1+wrap)%3],lookup[(2+wrap)%3]};
                        int wc[3]={rc[0],rc[1],rc[2]};
                        for(unsigned int swap = 0 ; swap < 2 ; swap++)
                        {
                            if(swap==1)
                            {
                                int el2=rc[1];
                                int el3=rc[2];
                                rc[2]=el2;
                                rc[1]=el3;
                            }
                            //change the signs of each element
                            for(unsigned int a = 0 ; a < 2 ; a++)
                            {
                                for(unsigned int b = 0 ; b < 2 ; b++)
                                {
                                    for(unsigned int c = 0 ; c < 2 ; c++)
                                    {
                                        wc[0]=rc[0]*pow(-1,a+1);
                                        wc[1]=rc[1]*pow(-1,b+1);
                                        wc[2]=rc[2]*pow(-1,c+1);
                                        int lookupvec[3]={wc[0]+mypos[0],wc[1]+mypos[1],wc[2]+mypos[2]};
                                        //std::cout << lookupvec[0] << "\t" << lookupvec[1] << "\t" << lookupvec[2] << std::endl;
                                        unsigned int check_lookup=0;
                                        for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                                        {
                                            if(lookupvec[xyz]<0 && config::pbc[xyz]==true)
                                            {
                                                lookupvec[xyz]=geom::dim[xyz]*geom::Nk[xyz]+lookupvec[xyz];
                                                check_lookup++;
                                            }
                                            else if(lookupvec[xyz]>=geom::dim[xyz]*geom::Nk[xyz] && config::pbc[xyz]==true)
                                            {
                                                lookupvec[xyz]=lookupvec[xyz]-geom::dim[xyz]*geom::Nk[xyz];
                                                check_lookup++;
                                            }
                                            else if(lookupvec[xyz]>=0 && lookupvec[xyz]<geom::dim[xyz]*geom::Nk[xyz])
                                            {
                                                check_lookup++;
                                            }
                                        }

                                        //if check_lookup==3 then the lookup for the atom exists
                                        if(check_lookup==3)//then we are array safe (i.e. should not seg fault) to look up the atom
                                        {
                                            //find the species of the neighbouring atom
                                            unsigned int spec=geom::coords(lookupvec[0],lookupvec[1],lookupvec[2],1);
                                            //check if we have assigned a J value to this interaction vector previously
                                            //we also need to check that we are looking at the right species
                                            if(check(lookupvec[0],lookupvec[1],lookupvec[2])==0 && spec==s1)
                                            {
                                                if((config::exchm>98 && lookupvec[0]!=mypos[0]) || config::exchm<99)
                                                {
                                                    unsigned int neigh=geom::coords(lookupvec[0],lookupvec[1],lookupvec[2],0);
                                                    //in this case we have not looked it up before



                                                    //This is the adjncy (neighbour lookup) in the CSR format
                                                    tadjncy.push_back(neigh);
                                                    //std::cout << tadjncy[adjncy_counter] << std::endl;
                                                    //std::cin.get();
                                                    tadjncy[adjncy_counter]=neigh;
                                                    adjncy_counter++;
                                                    //determine the diagonal
                                                    if(i>neigh)//then we are on the lower diagonal
                                                    {
                                                        int tmp=i,count=0;
                                                        while(neigh!=tmp){tmp--;count++;}
                                                        checkdiag(geom::nspins-count)++;
                                                    }
                                                    else if(i==neigh)
                                                    {
                                                        checkdiag[geom::nspins]++;
                                                    }
                                                    else if(i<neigh)
                                                    {
                                                        int tmp=neigh,count=0;
                                                        while(tmp!=i){tmp--;count++;}
                                                        checkdiag[geom::nspins+count]++;
                                                    }
                                                    for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                                                    {
                                                        for(unsigned int beta = 0 ; beta < 3 ; beta++)
                                                        {
                                                            tdata[alpha][beta].push_back(J(sl,spec,shell,alpha,beta)/(geom::ucm.GetMu(aiuc)*llg::muB));

                                                        }
                                                    }
                                                    if(outputJ)
                                                    {
                                                        opJ << i << "\t" << neigh << "\t" << tdata[0][0][adjncy_counter-1] << "\t[ " << mypos[0] << "," << mypos[1] << "," << mypos[2] << " ] -> [ " << lookupvec[0] << "," << lookupvec[1] << "," << lookupvec[2] << " ]" << std::endl;
                                                    }
                                                    check(lookupvec[0],lookupvec[1],lookupvec[2])=1;//this should mean that we no longer look for a neighbour here to avoid double addition of exchange
                                                }
                                            }
                                        }

                                    }
                                }
                            }
                        }

                    }
                }
            }
            xadj_counter++;
        }
        int diagcount=0;
        for(unsigned int i = 0 ; i < 2*geom::nspins-1 ; i++)
        {
            if(checkdiag[i]>0){diagcount++;}
        }
        diagnumdiag=diagcount;
        //std::cout <<"Number of diagonals detected\t" << diagcount << std::endl;
        xadj[geom::nspins]=adjncy_counter;
        if(outputJ)
        {
            opJ.close();
            if(opJ.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errWarning("Could not close J.dat for writing interaction matrix.");
            }
        }
        if(config::offdiag==false)
        {
            adjncy.resize(adjncy_counter);
            dataxx.resize(adjncy_counter);
            datayy.resize(adjncy_counter);
            datazz.resize(adjncy_counter);
            for(unsigned int i = 0 ; i < adjncy_counter ; i++)
            {
                adjncy[i]=tadjncy[i];
                dataxx[i]=tdata[0][0][i];
                datayy[i]=tdata[1][1][i];
                datazz[i]=tdata[2][2][i];
            }
            //clear the temporary arrays to free up memory
            tdata.clear();
            tadjncy.clear();
            FIXOUT(config::Info,"Size of xadj (CSR offsets):\t" << xadj.size() << std::endl);
            FIXOUT(config::Info,"Size of adjncy (CSR list):\t" << adjncy.size() << std::endl);
        }
    }
}
