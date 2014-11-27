// File: exch.cpp
// Author: Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 27 Nov 2014 10:48:54
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
    bool outputJ;
    void initExch(int argc,char *argv[])
    {

        std::string readMethod,readFile,method;
        libconfig::Config exchcfg;
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Exchange details***" << std::endl;
        try
        {
            config::cfg.readFile(argv[1]);
        }
        catch(const libconfig::FileIOException &fioex)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("I/O error while reading config file");
        }
        catch(const libconfig::ParseException &pex)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::cerr << ". Parse error at " << pex.getFile()  << ":" << pex.getLine() << "-" << pex.getError() << "***\n" << std::endl;
            exit(EXIT_FAILURE);
        }
        libconfig::Setting &setting = config::cfg.lookup("exchange");
        setting.lookupValue("exchmethod",method);
        FIXOUT(config::Info,"Method to read in exchange:" << method << std::endl);
        setting.lookupValue("exchinput",readMethod);
        setting.lookupValue("OutputExchange",outputJ);
        FIXOUT(config::Info,"Output of exchange matrix (mostly for visualization how diagonal it is):" << config::isTF(outputJ) << std::endl);
        if(readMethod=="thisfile")
        {
            FIXOUT(config::Info,"Reading exchange interactions from:" << "config file" << std::endl);
        }
        else if(readMethod=="extfile")
        {
            FIXOUT(config::Info,"Reading exchange interactions from:" << "external file" << std::endl);
            setting.lookupValue("exchfile",readFile);
            FIXOUT(config::Info,"File:" << readFile << std::endl);
            if(method=="permute")
            {
                try
                {
                    exchcfg.readFile(readFile.c_str());
                }
                catch(const libconfig::FileIOException &fioex)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("I/O error while reading exchange config file");
                }
                catch(const libconfig::ParseException &pex)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::cerr << ". Parse error at " << pex.getFile()  << ":" << pex.getLine() << "-" << pex.getError() << "***\n" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Method of reading exchange file not recognised.");
        }
        libconfig::Setting &GlobExch = exchcfg.lookup("exchange");
        GlobExch.lookupValue("MaxShells",max_shells);
        evs.resize(3);
        evs[0]=GlobExch["Scale"][0];
        evs[1]=GlobExch["Scale"][1];
        evs[2]=GlobExch["Scale"][2];
        FIXOUTVEC(config::Info,"Scaling factor for exchange vectors:",evs[0],evs[1],evs[2]);
        numint.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells);
        exchvec.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells,3);
        shell_list.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS());
        numint.IFill(0);exchvec.IFill(0);J.IFill(0);shell_list.IFill(0);
        if(geom::nspins>1)//then really we don't want any exchange anyway.
        {
            if(method=="permute")
            {
                //first read the exchange constants
                J.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells,3,3);
                for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
                {
                    for(unsigned int s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
                    {

                        config::printline(config::Info);
                        FIXOUT(config::Info,"Exchange interaction between species:" << s1 << " and " << s2 << std::endl);
                        std::stringstream sstr_int;
                        sstr_int << "exchange" << "_" << s1 << "_" << s2;
                        std::string str_int=sstr_int.str();
                        libconfig::Setting &exchset = exchcfg.lookup(str_int.c_str());
                        exchset.lookupValue("Num_Shells",shell_list(s1,s2));
                        FIXOUT(config::Info,"Reading exchange information for:" << shell_list(s1,s2) << " shells" << std::endl);
                        //Read the units of the exchange energy and scale appropriately to Joules
                        exchset.lookupValue("units",enerType);
                        FIXOUT(config::Info,"Units of exchange:" << enerType << std::endl);

                        for(unsigned int i = 0 ; i < shell_list(s1,s2) ; i++)
                        {
                            std::stringstream nisstr,evsstr;
                            std::string nistr,evstr;
                            nisstr << "NumInt" << i+1;
                            nistr=nisstr.str();
                            evsstr << "Shell" << i+1 << "Vec";
                            evstr=evsstr.str();
                            exchset.lookupValue(nistr.c_str(),numint(s1,s2,i));
                            config::Info << "Shell " << i+1;
                            FIXOUT(config::Info," number of interactions:" << numint(s1,s2,i) << std::endl);
                            for(unsigned int j = 0 ; j < 3 ; j++)
                            {
                                exchvec(s1,s2,i,j)=exchset[evstr.c_str()][j];
                            }
                            FIXOUTVEC(config::Info,"Vectors:",exchvec(s1,s2,i,0),exchvec(s1,s2,i,1),exchvec(s1,s2,i,2));
                            config::Info << std::endl;
                            for(unsigned int j = 0 ; j < 3 ; j++)
                            {

                                std::stringstream Jsstr;
                                Jsstr << "J" << i+1 << "_" << j+1;
                                std::string Jstr;
                                Jstr=Jsstr.str();
                                for(unsigned int k = 0 ; k < 3 ; k++)
                                {
                                    J(s1,s2,i,j,k) = exchset[Jstr.c_str()][k];
                                    if(enerType=="mRy")
                                    {
                                        J(s1,s2,i,j,k)*=2.179872172e-18; //now to milli
                                        J(s1,s2,i,j,k)*=1.0e-3;
                                    }
                                    else if(enerType=="eV")
                                    {
                                        J(s1,s2,i,j,k)*=1.602176565e-19;
                                    }
                                    else if(enerType=="J" || enerType=="Joules" || enerType=="joules")
                                    {
                                        //do nothing
                                    }
                                    else
                                    {
                                        error::errPreamble(__FILE__,__LINE__);
                                        error::errMessage("Units of exchange energy not recognised");
                                    }

                                }
                                FIXOUTVEC(config::Info,Jstr,J(s1,s2,i,j,0),J(s1,s2,i,j,1),J(s1,s2,i,j,2));
                            }

                            config::Info << std::endl;

                        }
                    }
                }
                checkdiag.resize((geom::nspins*2)-1);
                checkdiag.IFill(0);
                //then add the exchange constants to the matrix or to the interaction matrix
                if(config::exchm>0)//We calculate the exchange via a matrix multiplication (or part of it at least)
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
                                for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                                {
                                    //reference array
                                    int rc[3]={lookup[wrap%3],lookup[(1+wrap)%3],lookup[(2+wrap)%3]};
                                    int wc[3]={rc[0],rc[1],rc[2]};
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
                        //For debugging the CSR neighbour list
                        /*for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            std::cout << "Atom " << i << " has neighbours" << std::endl;
                            for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
                            {
                                std::cout << adjncy[j] << "\t" << dataxx[j] << "\t" << std::flush;
                            }
                            std::cout << std::endl;
                        }
                        std::cin.get();*/

                        if(config::exchm==1)
                        {//then convert CSR to DIA

                            //call the routine to convert the CSR to a DIA sparse matrix format
                            FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                            matconv::csr_to_dia_diag(geom::nspins,xadj,adjncy,diagoffset,diagnumdiag,checkdiag,dataxx,datayy,datazz);
                            SUCCESS(config::Info);
                            FIXOUT(config::Info,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                            config::openLogFile();
                            config::printline(config::Log);
                            FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                            config::Log << " [ ";
                            for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                            {
                                config::Log << diagoffset[i] << ",";
                            }
                            config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                            //for debugging the data array
                            /*for(int i = 0 ; i < diagnumdiag ; i++)
                              {
                              int os=diagoffset[i];
                              std::cout << "Offset " << i << " is " << os << std::endl;
                              for(int j = 0 ; j < geom::nspins ; j++)
                              {
                              std::cout << j << "\t" << i*geom::nspins+j << "\t" << datazz[i*geom::nspins+j] << std::endl;
                              }
                              std::cin.get();
                              }*/
                        }
                        #ifdef CUDA
                        if(config::exchm>98)
                        {
                            //call the routine to convert the CSR to a DIA sparse matrix format
                            FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                            matconv::csr_to_dia_diag(geom::nspins,xadj,adjncy,diagoffset,diagnumdiag,checkdiag,dataxx,datayy,datazz);
                            SUCCESS(config::Info);
                            FIXOUT(config::Info,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                            config::openLogFile();
                            config::printline(config::Log);
                            FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                            config::Log << " [ ";
                            for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                            {
                                config::Log << diagoffset[i] << ",";
                            }
                            config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                        }
                        if(config::offdiag && config::exchm>98)
                        {
                            offdiagoffset.resize(diagoffset.size());
                            for(unsigned int i = 0 ; i < offdiagoffset.size() ; i++)
                            {
                                offdiagoffset[i]=diagoffset[i];
                            }
                            FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                            matconv::csr_to_dia_offdiag(geom::nspins,xadj,adjncy,offdiagoffset,diagnumdiag,checkdiag,dataxy,dataxz,datayx,datayz,datazx,datazy);
                            SUCCESS(config::Info);
                            FIXOUT(config::Log,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                            config::openLogFile();
                            config::printline(config::Log);
                            FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                            config::Log << " [ ";
                            for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                            {
                                config::Log << diagoffset[i] << ",";
                            }
                            config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                            FIXOUT(config::Log,"Total number of non-zero diagonals for the anti-symmetric part of the exchange tensor:" << offdiagnumdiag << std::endl);
                            config::printline(config::Log);
                            FIXOUT(config::Log,"Outputting offset for off-diagonal (antisym) part of interaction matrix:" << std::endl);
                            config::Log << " [ ";
                            for(unsigned int i = 0 ; i < offdiagnumdiag-1 ; i++)
                            {
                                config::Log << offdiagoffset[i] << ",";
                            }
                            config::Log << offdiagoffset[offdiagnumdiag-1] << " ] " << std::endl;


                        }
                        #endif /*CUDA*/
                        if(config::offdiag && config::exchm==1)
                        {
                            offdiagoffset.resize(diagoffset.size());
                            for(unsigned int i = 0 ; i < offdiagoffset.size() ; i++)
                            {
                                offdiagoffset[i]=diagoffset[i];
                            }
                            FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                            matconv::csr_to_dia_offdiag(geom::nspins,xadj,adjncy,offdiagoffset,diagnumdiag,checkdiag,dataxy,dataxz,datayx,datayz,datazx,datazy);
                            SUCCESS(config::Info);
                            FIXOUT(config::Log,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                            config::openLogFile();
                            config::printline(config::Log);
                            FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                            config::Log << " [ ";
                            for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                            {
                                config::Log << diagoffset[i] << ",";
                            }
                            config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                            FIXOUT(config::Log,"Total number of non-zero diagonals for the anti-symmetric part of the exchange tensor:" << offdiagnumdiag << std::endl);
                            config::printline(config::Log);
                            FIXOUT(config::Log,"Outputting offset for off-diagonal (antisym) part of interaction matrix:" << std::endl);
                            config::Log << " [ ";
                            for(unsigned int i = 0 ; i < offdiagnumdiag-1 ; i++)
                            {
                                config::Log << offdiagoffset[i] << ",";
                            }
                            config::Log << offdiagoffset[offdiagnumdiag-1] << " ] " << std::endl;


                        }
                        if(config::exchm!=2)//then we are not using CSR
                        {
                            xadj.clear();
                            adjncy.clear();
                            checkdiag.clear();
                        }
                    }
                }
                else if(config::exchm==0)//add the exchange to the interaction matrix
                {
//check if we have a single layer of atoms in any dimension
                    bool checkmonolayer[3]={false,false,false};
                    for(unsigned int xyz = 0 ; xyz < 3; xyz++)
                    {
                        if(geom::dim[xyz]*geom::Nk[xyz] < 2)
                        {
                            checkmonolayer[xyz]=true;
                        }
                    }
                    FIXOUTVEC(config::Info,"Monolayer check:",config::isTF(checkmonolayer[0]),config::isTF(checkmonolayer[1]),config::isTF(checkmonolayer[2]));

                    //we are going to write the exchange information to the log file. Make sure it if open.
                    config::openLogFile();
                    for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
                    {
                        config::Log << "Exchange parameters acting on species " << s1 << std::endl;
                        for(unsigned int s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
                        {
                            config::Log << "Interaction with species " << s2 << std::endl;
                            Array3D<unsigned int> check;
                            check.resize(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
                            check.IFill(0);
                            //This section of code takes the kvec interactions and
                            //adds the appropriate Jij to the appropriate interaction matrix
                            for(unsigned int i = 0 ; i < shell_list(s1,s2) ; i++)
                            {
                                config::Log << "Shell " << shell_list(s1,s2) << " of " << numint(s1,s2,i) << std::endl;
                                unsigned int counter=0;
                                int lc[3]={0,0,0};
                                lc[0]=static_cast<unsigned int>(exchvec(s1,s2,i,0)*static_cast<double>(geom::Nk[0])*evs[0]+0.5);
                                lc[1]=static_cast<unsigned int>(exchvec(s1,s2,i,1)*static_cast<double>(geom::Nk[1])*evs[1]+0.5);
                                lc[2]=static_cast<unsigned int>(exchvec(s1,s2,i,2)*static_cast<double>(geom::Nk[2])*evs[2]+0.5);
                                for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                                {
                                    //reference array
                                    int rc[3]={lc[wrap%3],lc[(1+wrap)%3],lc[(2+wrap)%3]};
                                    //work array
                                    int wc[3]={rc[0],rc[1],rc[2]};
                                    if((abs(wc[0]>0) && checkmonolayer[0]==true) || (abs(wc[1]>0) && checkmonolayer[1]==true) || (abs(wc[2]>0) && checkmonolayer[2]==true) )
                                    {
                                        //then do nothing, we don't want to add anything in this direction
                                    }
                                    else
                                    {
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
                                                    //This array is purely for outputting the exchange information to the log file
                                                    int owc[3]={wc[0],wc[1],wc[2]};
                                                    //check the boundaries for each component
                                                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                                                    {
                                                        if(wc[xyz]<0)
                                                        {
                                                            wc[xyz]=geom::zpdim[xyz]*geom::Nk[xyz]+wc[xyz];
                                                        }
                                                    }//end of xyz loop
                                                    if(check(wc[0],wc[1],wc[2])==0)
                                                    {
                                                        config::Log << "Interaction Vector:  [" << owc[0] << "," << owc[1] << "," << owc[2] << "]\t -> [" << wc[0] << "," << wc[1] << "," << wc[2] << "]" << std::endl;
                                                        //add the diagonal components (alpha=beta)
                                                        for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                                                        {
                                                            intmat::Nrab(s1,s2,alpha,alpha,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,alpha,alpha)/(geom::ucm.GetMuBase(s1)*llg::muB));
                                                        }
                                                        //do the DM (off-diagonals by hand)
                                                        //The format of the file that is read in is in Jxx. We want in our interaction
                                                        //matrix the DM vectors.
                                                        // Nxy = 1/2(Jyx-Jxy)
                                                        intmat::Nrab(s1,s2,0,1,wc[0],wc[1],wc[2])[0]+=(0.5*(J(s1,s2,i,1,0)-J(s1,s2,i,0,1)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                        // Nxz = 1/2(Jxz-Jzx)
                                                        intmat::Nrab(s1,s2,0,2,wc[0],wc[1],wc[2])[0]+=(0.5*(J(s1,s2,i,0,2)-J(s1,s2,i,2,0)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                        // Nyx = 1/2(Jxy-Jyx)
                                                        intmat::Nrab(s1,s2,1,0,wc[0],wc[1],wc[2])[0]+=(0.5*(J(s1,s2,i,0,1)-J(s1,s2,i,1,0)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                        // Nyz = 1/2(Jzy-Jyz)
                                                        intmat::Nrab(s1,s2,1,2,wc[0],wc[1],wc[2])[0]+=(0.5*(J(s1,s2,i,2,1)-J(s1,s2,i,1,2)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                        // Nzx = 1/2(Jzx - Jxz)
                                                        intmat::Nrab(s1,s2,2,0,wc[0],wc[1],wc[2])[0]+=(0.5*(J(s1,s2,i,2,0)-J(s1,s2,i,0,2)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                        // Nzy = 1/2(Jyz-Jzy)
                                                        intmat::Nrab(s1,s2,2,1,wc[0],wc[1],wc[2])[0]+=(0.5*(J(s1,s2,i,1,2)-J(s1,s2,i,2,1)))/(geom::ucm.GetMuBase(s1)*llg::muB);

                                                        config::Log << "[ " << J(s1,s2,i,0,0) << " , " << J(s1,s2,i,0,1) << " , " << J(s1,s2,i,0,2) << " ]" << std::endl;
                                                        config::Log << "[ " << J(s1,s2,i,1,0) << " , " << J(s1,s2,i,1,1) << " , " << J(s1,s2,i,1,2) << " ]\t (Joules)" << std::endl;
                                                        config::Log << "[ " << J(s1,s2,i,2,0) << " , " << J(s1,s2,i,2,1) << " , " << J(s1,s2,i,2,2) << " ]" << std::endl;
                                                        config::Log << std::endl;
                                                        config::Log << "[ " << J(s1,s2,i,0,0)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s2,i,0,1)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s2,i,0,2)/(geom::ucm.GetMuBase(s1)*llg::muB) << " ]" << std::endl;
                                                        config::Log << "[ " << J(s1,s2,i,1,0)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s2,i,1,1)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s2,i,1,2)/(geom::ucm.GetMuBase(s1)*llg::muB) << " ]\t (Tesla)" << std::endl;
                                                        config::Log << "[ " << J(s1,s2,i,2,0)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s2,i,2,1)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s2,i,2,2)/(geom::ucm.GetMuBase(s1)*llg::muB) << " ]" << std::endl;
                                                        config::Log << std::endl;
                                                        check(wc[0],wc[1],wc[2])=1;
                                                        counter++;

                                                    }//end of check if statement
                                                }//end of c loop
                                            }//end of b loop
                                        }//end of a loop
                                    }//end of monolayer if statement
                                }//end of wrap loop
                                if(counter!=numint(s1,s2,i))
                                {
                                    error::errPreamble(__FILE__,__LINE__);
                                    std::stringstream sstr;
                                    sstr << "The number of interactions between species " << s1 << " and " << s2 << " in shell " << shell_list(s1,s2) << " should be " << numint(s1,s2,i) << ", instead the number counted was " << counter;
                                    std::string str=sstr.str();
                                    error::errMessage(str.c_str());
                                }//end of counter check if statement
                                counter=0;
                                config::printline(config::Log);
                            }//end of shell list loop
                        }//end of s2 loop
                    }//end of s1 loop
                    /*                for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
                                      {
                                      for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
                                      {
                                      for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                                      {
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(0,0,2,2,i,j,k)[0] << std::endl;
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(0,1,2,2,i,j,k)[0] << std::endl;
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(1,0,2,2,i,j,k)[0] << std::endl;
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(1,1,2,2,i,j,k)[0] << std::endl;
                                      std::cin.get();
                                      }
                                      }
                                      }
                     */
                    /*                for(unsigned int i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
                                      {
                                      for(unsigned int j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
                                      {
                                      for(unsigned int k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                                      {
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(0,0,2,2,i,j,k)[0] << std::endl;
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(0,1,2,2,i,j,k)[0] << std::endl;
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(1,0,2,2,i,j,k)[0] << std::endl;
                                      std::cout << i << "\t" << j << "\t" << k << "\t" << intmat::Nrab(1,1,2,2,i,j,k)[0] << std::endl;
                                      std::cin.get();
                                      }
                                      }
                                      }
                     */
                }


                if(config::exchm>98)
                {
                    bool checkmonolayer[3]={false,false,false};
                    for(unsigned int xyz = 0 ; xyz < 3; xyz++)
                    {
                        if(geom::dim[xyz]*geom::Nk[xyz] < 2)
                        {
                            checkmonolayer[xyz]=true;
                        }
                    }
                    FIXOUTVEC(config::Info,"Monolayer check:",config::isTF(checkmonolayer[0]),config::isTF(checkmonolayer[1]),config::isTF(checkmonolayer[2]));

                    //we are going to write the exchange information to the log file. Make sure it if open.
                    config::openLogFile();
                    Array2D<unsigned int> check;
                    check.resize(geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);
                    for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)//here geom::ucm.GetNMS() should be equal to the number of planes (dim[0]*Nk[0])
                    {
                        //for checking eqch plane we need to reset the check array
                        check.IFill(0);
                        //This section of code takes the kvec interactions and
                        //adds the appropriate Jij to the appropriate interaction matrix
                        for(unsigned int i = 0 ; i < shell_list(s1,s1) ; i++)
                        {
                            config::Log << "Shell " << shell_list(s1,s1) << " of " << numint(s1,s1,i) << std::endl;
                            unsigned int counter=0;
                            int lc[3]={0,0,0};
                            int plane=s1;
                            if(s1 > geom::dim[0]*geom::Nk[0])
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("The number of species should be equal to the number of planes in the x-direcion to use the hybrid method.");
                            }

                            lc[0]=static_cast<unsigned int>(exchvec(s1,s1,i,0)*static_cast<double>(geom::Nk[0])*evs[0]+0.5);;
                            lc[1]=static_cast<unsigned int>(exchvec(s1,s1,i,1)*static_cast<double>(geom::Nk[1])*evs[1]+0.5);
                            lc[2]=static_cast<unsigned int>(exchvec(s1,s1,i,2)*static_cast<double>(geom::Nk[2])*evs[2]+0.5);
                            for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                            {
                                //reference array
                                int rc[3]={lc[(0+wrap)%3],lc[(1+wrap)%3],lc[(2+wrap)%3]};

                                //work array
                                int wc[3]={rc[0],rc[1],rc[2]};

                                if( (abs(wc[0]>0) && checkmonolayer[0]==true) || (abs(wc[1]>0) && checkmonolayer[1]==true) || (abs(wc[2]>0) && checkmonolayer[2]==true))
                                {
                                    //then do nothing, we don't want to add anything in this direction
                                }
                                else
                                {

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

                                                //This array is purely for outputting the exchange information to the log file
                                                int owc[3]={wc[0],wc[1],wc[2]};
                                                //check the boundaries for each component
                                                for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                                                {
                                                    if(wc[xyz]<0)
                                                    {
                                                        wc[xyz]=geom::zpdim[xyz]*geom::Nk[xyz]+wc[xyz];
                                                    }
                                                }//end of xyz loop
                                                unsigned int sum=abs(wc[0])+abs(wc[1])+abs(wc[2]);
                                                if(check(wc[1],wc[2])==0 && wc[0]==0 && sum>0)
                                                {
                                                    config::Log << "Plane " << plane << " Interaction Vector:  [" << owc[1] << "," << owc[2] << "]\t -> [" << wc[0] << "," << wc[1] << "," << wc[2] << "]" << std::endl;
                                                    //add the diagonal components (alpha=beta)
                                                    for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                                                    {
                                                        intmat::hNrab(alpha,alpha,plane,wc[1],wc[2])[0]+=(J(s1,s1,i,alpha,alpha)/(geom::ucm.GetMuBase(s1)*llg::muB));
//                                                        if(alpha==0)
//                                                        {
//                                                            std::cout << "Plane Interaction Vector:  [" << plane << "," << owc[1] << "," << owc[2] << "]\t -> [" << wc[0] << "," << wc[1] << "," << wc[2] << "]" << std::endl;
//                                                            std::cout << "Adding to intmat:\t" << intmat::hNrab(alpha,alpha,plane,wc[1],wc[2])[0] << std::endl;
//                                                        }
                                                    }
                                                    //do the DM (off-diagonals by hand)
                                                    //The format of the file that is read in is in Jxx. We want in our interaction
                                                    //matrix the DM vectors.
                                                    // Nxy = 1/2(Jyx-Jxy)
                                                    intmat::hNrab(0,1,plane,wc[1],wc[2])[0]+=(0.5*(J(s1,s1,i,1,0)-J(s1,s1,i,0,1)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                    // Nxz = 1/2(Jxz-Jzx)
                                                    intmat::hNrab(0,2,plane,wc[1],wc[2])[0]+=(0.5*(J(s1,s1,i,0,2)-J(s1,s1,i,2,0)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                    // Nyx = 1/2(Jxy-Jyx)
                                                    intmat::hNrab(1,0,plane,wc[1],wc[2])[0]+=(0.5*(J(s1,s1,i,0,1)-J(s1,s1,i,1,0)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                    // Nyz = 1/2(Jzy-Jyz)
                                                    intmat::hNrab(1,2,plane,wc[1],wc[2])[0]+=(0.5*(J(s1,s1,i,2,1)-J(s1,s1,i,1,2)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                    // Nzx = 1/2(Jzx - Jxz)
                                                    intmat::hNrab(2,0,plane,wc[1],wc[2])[0]+=(0.5*(J(s1,s1,i,2,0)-J(s1,s1,i,0,2)))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                                    // Nzy = 1/2(Jyz-Jzy)
                                                    intmat::hNrab(2,1,plane,wc[1],wc[2])[0]+=(0.5*(J(s1,s1,i,1,2)-J(s1,s1,i,2,1)))/(geom::ucm.GetMuBase(s1)*llg::muB);

                                                    config::Log << "[ " << J(s1,s1,i,0,0) << " , " << J(s1,s1,i,0,1) << " , " << J(s1,s1,i,0,2) << " ]" << std::endl;
                                                    config::Log << "[ " << J(s1,s1,i,1,0) << " , " << J(s1,s1,i,1,1) << " , " << J(s1,s1,i,1,2) << " ]\t (Joules)" << std::endl;
                                                    config::Log << "[ " << J(s1,s1,i,2,0) << " , " << J(s1,s1,i,2,1) << " , " << J(s1,s1,i,2,2) << " ]" << std::endl;
                                                    config::Log << std::endl;
                                                    config::Log << "[ " << J(s1,s1,i,0,0)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s1,i,0,1)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s1,i,0,2)/(geom::ucm.GetMuBase(s1)*llg::muB) << " ]" << std::endl;
                                                    config::Log << "[ " << J(s1,s1,i,1,0)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s1,i,1,1)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s1,i,1,2)/(geom::ucm.GetMuBase(s1)*llg::muB) << " ]\t (Tesla)" << std::endl;
                                                    config::Log << "[ " << J(s1,s1,i,2,0)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s1,i,2,1)/(geom::ucm.GetMuBase(s1)*llg::muB) << " , " << J(s1,s1,i,2,2)/(geom::ucm.GetMuBase(s1)*llg::muB) << " ]" << std::endl;
                                                    config::Log << std::endl;
                                                    check(wc[1],wc[2])=1;
                                                    counter++;

                                                }//end of check if statement
                                            }//end of c loop
                                        }//end of b loop
                                    }//end of a loop
                                }//end of wrap loop
                            }
                                counter=0;
                                config::printline(config::Log);
                            }//end of shell list loop
                        }//end of s1 loop


                    }//end of else if(config::exchm>98)

                }//end of if(method=="permute") statement
            else if(method=="direct")
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("The use of \"direct\" method of reading in the exchange is currently not written for a general N lattice code and needs to be fixed");

                std::ifstream ifs;
                ifs.open(readFile.c_str());
                if(!ifs.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open exchange file for reading");
                }
                unsigned int noint = 0;
                ifs >> noint;
                FIXOUT(config::Info,"Number of interactions to be read in:" << noint << std::endl);
                Array3D<unsigned int> check(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);//This array is used to check if we already have the Jij for this interaction
                check.IFill(0);
                unsigned int counter=0;

                int dump;
                ifs>>dump;
                std::ofstream map("map.dat");
                if(!map.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open file for outputting exchange map");
                }
                double dcut=5.01;
                for(unsigned int i = 0 ; i < noint ; i++)
                {
                    int oc[3]={0,0,0},c[3]={0,0,0};
                    double J[3][3];
                    ifs>>dump;
                    ifs>>dump;
                    ifs>>dump;
                    double dist=0.0;
                    for(unsigned int rc = 0 ; rc < 3 ; rc++)
                    {
                        ifs >> c[rc];
                        oc[rc]=c[rc];

                        dist+=(static_cast<double>(c[rc])*static_cast<double>(c[rc]));
                        //c[rc]*=geom::Nk[rc];
                        //check boundaries
                        if(c[rc]<0)
                        {
                            c[rc]=geom::zpdim[rc]*geom::Nk[rc]+c[rc];
                        }
                    }
                    dist=sqrt(dist);

                    //std::cout << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << 1 << std::endl;
                    for(unsigned int j1 = 0 ; j1 < 3 ; j1++)
                    {
                        for(unsigned int j2 = 0 ; j2 < 3 ; j2++)
                        {
                            ifs >> J[j1][j2];
                        }
                    }
                    if(c[2]==0)
                    {
                        if(sqrt(oc[0]*oc[0]*0.5*0.5+oc[1]*sqrt(2.0)*oc[1]*sqrt(2.0)/16) <dcut)
                        {
                            map << oc[0] << "\t" << oc[1];
                            for(unsigned int j1 = 0 ; j1 < 3 ; j1++)
                            {
                                for(unsigned int j2 = 0 ; j2 < 3 ; j2++)
                                {
                                    map << "\t" << J[j1][j2];
                                }
                            }
                            map << std::endl;
                        }
                    }

                    if(check(c[0],c[1],c[2])==0)//then we do not already have an interaction there
                    {
                        check(c[0],c[1],c[2])=1;
                        counter++;
                        if(geom::coords(c[0],c[1],c[2],0)>-2)
                        {
                            if(sqrt(oc[0]*oc[0]*0.5*0.5+oc[1]*sqrt(2.0)*oc[1]*sqrt(2.0)/16) <dcut)
                            {
                                //intmat::Nrxx(c[0],c[1],c[2])+=(J[0][0]/(mat::muB*mat::mu));
                                //intmat::Nryy(c[0],c[1],c[2])+=(J[1][1]/(mat::muB*mat::mu));
                                //intmat::Nrzz(c[0],c[1],c[2])+=(J[2][2]/(mat::muB*mat::mu));
                                //The format of the file that is read in is in Jxx. We want in our interaction
                                //matrix the DM vectors.
                                // Nxy = 1/2(Jyx-Jxy)
                                //intmat::Nrxy(c[0],c[1],c[2])+=((0.5*(J[1][0]-J[0][1]))/(mat::muB*mat::mu));
                                // Nxz = 1/2(Jxz-Jzx)
                                //intmat::Nrxz(c[0],c[1],c[2])+=((0.5*(J[0][2]-J[2][0]))/(mat::muB*mat::mu));
                                // Nyx = 1/2(Jxy-Jyx)
                                //intmat::Nryx(c[0],c[1],c[2])+=((0.5*(J[0][1]-J[1][0]))/(mat::muB*mat::mu));
                                // Nyz = 1/2(Jzy-Jyz)
                                //intmat::Nryz(c[0],c[1],c[2])+=((0.5*(J[2][1]-J[1][2]))/(mat::muB*mat::mu));
                                // Nzx = 1/2(Jzx - Jxz)
                                //intmat::Nrzx(c[0],c[1],c[2])+=((0.5*(J[2][0]-J[0][2]))/(mat::muB*mat::mu));
                                // Nzy = 1/2(Jyz-Jzy)
                                //intmat::Nrzy(c[0],c[1],c[2])+=((0.5*(J[1][2]-J[2][1]))/(mat::muB*mat::mu));
                            }
                            //intmat::Nrzz(c[0],c[1],c[2])+=((J[2][2]+2.0*2.0*1.6e-19*1e-3)/(mat::muB*mat::mu));
                            //                std::cout << "Interaction: " << i << "\nJij:\n" << intmat::Nrxx(c[0],c[1],c[2]) << "\t" << intmat::Nrxy(c[0],c[1],c[2]) << "\t" << intmat::Nrxz(c[0],c[1],c[2]) << std::endl;
                            //                std::cout <<  intmat::Nryx(c[0],c[1],c[2]) << "\t" << intmat::Nryy(c[0],c[1],c[2]) << "\t" << intmat::Nryz(c[0],c[1],c[2]) << std::endl;
                            //                std::cout <<  intmat::Nrzx(c[0],c[1],c[2]) << "\t" << intmat::Nrzy(c[0],c[1],c[2]) << "\t" << intmat::Nrzz(c[0],c[1],c[2]) << std::endl;
                        }
                        else
                        {
                            std::cout << oc[0] << "\t" << oc[1] << "\t" << oc[2] << "\t" <<  c[0] << "\t" << c[1] << "\t" << c[2] << std::endl;
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("You are trying to add an interaction to an empty mesh point.");
                        }
                    }
                    else
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("That interaction has already been read");
                    }

                }
                if(counter!=noint)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Incorrect number of interactions found");
                }
                else
                {
                    FIXOUT(config::Info,"Read in exchange data for:" << counter << " interactions" << std::endl);
                }
                check.clear();
                map.close();
                if(map.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errWarning("Could not close interaction map file.");
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
