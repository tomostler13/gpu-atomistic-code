// File: exch_unitcell.cpp
// Author: Tom Ostler
// Created: 05 Dec 2014
// Last-modified: 13 Sep 2016 20:10:39
// This routine determines the exchange matrix for the unitcell method
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
    void unitcell(int argc,char *argv[])
    {
        get_exch_unitcell(argc,argv);
        checkdiag.resize((geom::nspins*2)-1);
        checkdiag.IFill(0);
        if(config::exchm>0)//We calculate the exchange via a matrix multiplication (or part of it at least)
        {

            if(config::exchm==2)
            {
                exchunitcell();
                //output the sparse matrix
                if(oem)
                {
                    if(config::offdiag==false)
                    {

                        std::ofstream opem("csr_exch_mat.dat");
                        if(!opem.is_open())
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("Could not open file for outputting dia exchange matrix.");
                        }
                        opem << geom::nspins << std::endl;
                        opem << adjncy.size() << std::endl;
                        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            opem << xadj[i] << "\t";
                        }
                        opem << xadj[geom::nspins] << std::endl;
                        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
                            {
                                opem << adjncy[j] << "\t" << dataxx[j] << "\t" << datayy[j] << "\t" << datazz[j] << "\t";
                            }
                        }
                        opem.close();
                        if(opem.is_open())
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errWarning("Could not close the file for outputting the CSR exchange matrix.");
                        }

                    }
                    else
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("FUTURE IMPLEMENTATION: At the moment the outputting of the off-diagonal CSR format is not possible.");
                    }

                }
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
                    if(oem)
                    {
                        if(config::offdiag)
                        {
                            std::ofstream opem("dia_exch_mat.dat");
                            if(!opem.is_open())
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("Could not open file for outputting dia exchange matrix.");
                            }
                            opem << geom::nspins << std::endl;
                            opem << diagnumdiag << std::endl;
                            for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                            {
                                opem << diagoffset[i] << "\t";
                            }
                            opem << diagoffset[diagnumdiag-1];
                            for(unsigned int i = 0 ; i < diagnumdiag ; i++)
                            {
                                for(unsigned int j = 0 ; j < geom::nspins ; j++)
                                {
                                    opem << dataxx[i*geom::nspins+j] << "\t" << datayy[i*geom::nspins+j] << "\t" << datazz[i*geom::nspins+j] << "\t";
                                }
                            }
                            opem.close();
                            if(opem.is_open())
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errWarning("Could not close the file for outputting the DIA exchange matrix.");
                            }
                        }
                        else
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("FUTURE IMPLEMENTATION: At the moment outputting of the off-diagonal in DIA format is not possible.");
                        }

                    }

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
                    if(oem)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("FUTURE IMPLEMENTATION: At the moment it is not possible to output the exchange matrix of the hybrid format.");
                    }


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
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The use of the FFT method with unitcell method for the exchange is not yet implemented.");
        }

        if(config::exchm>98)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The use of the hybrid method when reading the exchange using the unitcell method is not yet implemented.");
        }//end of else if(config::exchm>98)
    }

    void exchunitcell()
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
            opJ << "#No. of spin i - No. of spin j - Jij [T] - Coords (i) - Coords(j) - Distance" << std::endl;
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

        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            xadj[xadj_counter]=adjncy_counter;
            //get my atom in the unit cell
            unsigned int aiuc=geom::lu(i,4);
            //get my unit cell coords
            int ucc[3]={geom::lu(i,5),geom::lu(i,6),geom::lu(i,7)};
            //loop over the exchange
            for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
            {
                for(unsigned int shell = 0 ; shell < shell_list(aiuc,s1) ; shell++)
                {

                    int lookup[3]={0,0,0};
                    if(exchvec(aiuc,s1,shell,0)<0)
                    {
                        lookup[0]=static_cast<int>(exchvec(aiuc,s1,shell,0)-0.5);
                    }
                    else
                    {
                        lookup[0]=static_cast<int>(exchvec(aiuc,s1,shell,0)+0.5);
                    }
                    if(exchvec(aiuc,s1,shell,1)<0)
                    {
                        lookup[1]=static_cast<int>(exchvec(aiuc,s1,shell,1)-0.5);
                    }
                    else
                    {
                        lookup[1]=static_cast<int>(exchvec(aiuc,s1,shell,1)+0.5);
                    }
                    if(exchvec(aiuc,s1,shell,2)<0)
                    {
                        lookup[2]=static_cast<int>(exchvec(aiuc,s1,shell,2)-0.5);
                    }
                    else
                    {
                        lookup[2]=static_cast<int>(exchvec(aiuc,s1,shell,2)+0.5);
                    }
                    //std::cout << "Unit Cell Lookup (before) " << lookup[0] << "\t" << lookup[1] << "\t" << lookup[2] << std::endl;
                    //std::cout << "Input: " << exchvec(aiuc,s1,shell,0) << "\t"<< exchvec(aiuc,s1,shell,1) << "\t" << exchvec(aiuc,s1,shell,2) << std::endl;
                    int lookupvec[3]={lookup[0]+ucc[0],lookup[1]+ucc[1],lookup[2]+ucc[2]};
                    int cluv[3]={lookupvec[0],lookupvec[1],lookupvec[2]};
                    unsigned int check_lookup=0;
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        if(lookupvec[xyz]<0 && config::pbc[xyz]==true)
                        {
                            lookupvec[xyz]=geom::dim[xyz]+lookupvec[xyz];
                            check_lookup++;
                        }
                        else if(lookupvec[xyz]>=geom::dim[xyz] && config::pbc[xyz]==true)
                        {
                            lookupvec[xyz]=lookupvec[xyz]-geom::dim[xyz];
                            check_lookup++;
                        }
                        else if(lookupvec[xyz]>=0 && lookupvec[xyz]<geom::dim[xyz])//looks redundant but takes into account non-pbc's
                        {
                            check_lookup++;
                        }
                    }
                    //std::cout << "Unit Cell Lookup " << lookupvec[0] << "\t" << lookupvec[1] << "\t" << lookupvec[2] << std::endl;
                    //std::cin.get();

                    //check the neighbouring spin is not greater than the interaction range
                    const double ci[3]={geom::rx[i],geom::ry[i],geom::rz[i]};
                    double cj[3]={0.0,0.0,0.0};
                    double xred[3]={geom::ucm.GetXred(s1,0),geom::ucm.GetXred(s1,1),geom::ucm.GetXred(s1,2)};
                    cj[0]=(xred[0]+static_cast<double>(cluv[0]))*geom::rprim(0,0)+(xred[1]+static_cast<double>(cluv[1]))*geom::rprim(1,0)+(xred[2]+static_cast<double>(cluv[2]))*geom::rprim(2,0);
                    cj[1]=(xred[0]+static_cast<double>(cluv[0]))*geom::rprim(0,1)+(xred[1]+static_cast<double>(cluv[1]))*geom::rprim(1,1)+(xred[2]+static_cast<double>(cluv[2]))*geom::rprim(2,1);
                    cj[2]=(xred[0]+static_cast<double>(cluv[0]))*geom::rprim(0,2)+(xred[1]+static_cast<double>(cluv[1]))*geom::rprim(1,2)+(xred[2]+static_cast<double>(cluv[2]))*geom::rprim(2,2);
                    const double dist=sqrt((ci[0]-cj[0])*(ci[0]-cj[0])+(ci[1]-cj[1])*(ci[1]-cj[1])+(ci[2]-cj[2])*(ci[2]-cj[2]));
                    //std::cout << dist << "\t" << rcut << std::endl;
                    if(cutexch)
                    {
                        if(dist>rcut/1e-10)
                        {
                            check_lookup=0;
                        }
                    }

                    //std::cin.get();
                    if(check_lookup==3)
                    {
                        unsigned int neighid=geom::atnolu(lookupvec[0],lookupvec[1],lookupvec[2],s1);
                        tadjncy.push_back(neighid);
                        tadjncy[adjncy_counter]=neighid;
                        adjncy_counter++;
                        if(i>neighid)//then we are on the lower diagonal
                        {
                            int tmp=i,count=0;
                            while(neighid!=tmp){tmp--;count++;}
                            checkdiag(geom::nspins-count)++;
                        }
                        else if(i==neighid)
                        {
                            checkdiag[geom::nspins]++;
                        }
                        else if(i<neighid)
                        {
                            int tmp=neighid,count=0;
                            while(tmp!=i){tmp--;count++;}
                            checkdiag[geom::nspins+count]++;
                        }
                        for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                        {
                            for(unsigned int beta = 0 ; beta < 3 ; beta++)
                            {
                                tdata[alpha][beta].push_back(J(aiuc,s1,shell,alpha,beta)/(geom::ucm.GetMu(aiuc)*llg::muB));

                            }
                        }
                        if(outputJ)
                        {
                            opJ << i << "\t" << aiuc << "\t" << neighid << "\t" << s1 << "\t" << tdata[0][0][adjncy_counter-1] << "\t[ " << ci[0] << "," << ci[1] << "," << ci[2] << " ] -> [ " << cj[0] << "," << cj[1] << "," << cj[2] << " ]" << "\t" << dist << std::endl;
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
        else if(config::offdiag==true)
        {
            adjncy.resize(adjncy_counter);
            dataxx.resize(adjncy_counter);
            dataxy.resize(adjncy_counter);
            dataxz.resize(adjncy_counter);
            datayx.resize(adjncy_counter);
            datayy.resize(adjncy_counter);
            datayz.resize(adjncy_counter);
            datazx.resize(adjncy_counter);
            datazy.resize(adjncy_counter);
            datazz.resize(adjncy_counter);
            for(unsigned int i = 0 ; i < adjncy_counter ; i++)
            {
                adjncy[i]=tadjncy[i];
                dataxx[i]=tdata[0][0][i];
                dataxy[i]=tdata[0][1][i];
                dataxz[i]=tdata[0][2][i];
                datayx[i]=tdata[1][0][i];
                datayy[i]=tdata[1][1][i];
                datayz[i]=tdata[1][2][i];
                datazx[i]=tdata[2][0][i];
                datazy[i]=tdata[2][1][i];
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
