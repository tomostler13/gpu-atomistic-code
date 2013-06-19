// File: exch.cpp
// Author: Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 18 Jun 2013 12:44:21
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/mat.h"
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
    unsigned int num_shells = 0;
    Array<unsigned int> numint;
    Array2D<double> exchvec;
    Array2D<unsigned int> kvec;
    Array3D<double> J;
    Array<unsigned int> xadj,adjncy;
    Array<double> Jxx,Jyy,Jzz;
    std::string enerType;
    bool pbc[3];
    void initExch(int argc,char *argv[])
    {
        std::string readMethod,readFile,method;
        libconfig::Config exchcfg;
        config::printline(config::Info);
/*        std::ofstream intmap("interaction_map.dat");
        if(!(intmap.is_open()))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open interaction map file");
        }*/
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
        FIXOUT(config::Info,"Exchange read in:" << method << std::endl);
        setting.lookupValue("exchinput",readMethod);
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
        if(config::useintmat)
        {
            if(method=="permute")
            {
                libconfig::Setting &exchset = exchcfg.lookup("exchange");
                exchset.lookupValue("Num_Shells",num_shells);
                assert(num_shells>0);
                FIXOUT(config::Info,"Reading exchange information for:" << num_shells << " shells" << std::endl);
                numint.resize(num_shells);
                exchvec.resize(num_shells,3);
                kvec.resize(num_shells,3);
                J.resize(num_shells,3,3);
                //Read the units of the exchange energy and scale appropriately to Joules
                exchset.lookupValue("units",enerType);
                FIXOUT(config::Info,"Units of exchange:" << enerType << std::endl);

                for(unsigned int i = 0 ; i < num_shells ; i++)
                {
                    std::stringstream nisstr,evsstr;
                    std::string nistr,evstr;
                    nisstr << "NumInt" << i+1;
                    nistr=nisstr.str();
                    evsstr << "Shell" << i+1 << "Vec";
                    evstr=evsstr.str();
                    exchset.lookupValue(nistr.c_str(),numint[i]);
                    config::printline(config::Info);
                    config::Info << "Shell " << i+1;
                    FIXOUT(config::Info," number of interactions:" << numint[i] << std::endl);
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {
                        exchvec(i,j)=exchset[evstr.c_str()][j];
                        kvec(i,j)=int(exchvec(i,j)*geom::Nk[j]);
                    }
                    FIXOUTVEC(config::Info,"Vectors:",exchvec(i,0),exchvec(i,1),exchvec(i,2));
                    FIXOUTVEC(config::Info,"On K-mesh:",kvec(i,0),kvec(i,1),kvec(i,2));
                    FIXOUT(config::Info,"Distance in k-points:" << sqrt(double(kvec(i,0)*kvec(i,0)+kvec(i,1)*kvec(i,1)+kvec(i,2)*kvec(i,2))));
                    config::Info << std::endl;
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {

                        std::stringstream Jsstr;
                        Jsstr << "J" << i+1 << "_" << j+1;
                        std::string Jstr;
                        Jstr=Jsstr.str();
                        for(unsigned int k = 0 ; k < 3 ; k++)
                        {
                            J(i,j,k) = exchset[Jstr.c_str()][k];
                            if(enerType=="mRy")
                            {
                                J(i,j,k)*=2.179872172e-18;
                                //now to milli
                                J(i,j,k)*=1.0e-3;
                            }
                            else if(enerType=="eV")
                            {
                                J(i,j,k)*=1.602176565e-19;
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
                        FIXOUTVEC(config::Info,Jstr,J(i,j,0),J(i,j,1),J(i,j,2));
                    }

                    config::Info << std::endl;

                }

                //This section of code takes the kvec interactions and
                //adds the appropriate Jij to the appropriate interaction matrix

                Array3D<unsigned int> check(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);//This array is used to check if we already have the Jij for this interaction
                check.IFill(0);
                for(unsigned int i = 0 ; i < num_shells ; i++)
                {
                    unsigned int counter=0;
                    //store the original vector
                    unsigned int lc[3]={kvec(i,0),kvec(i,1),kvec(i,2)};
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {

                        if(geom::zpcheck==false && abs(lc[xyz])>geom::dim[xyz]*geom::Nk[xyz]/2)
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("Interactions are out of range");
                        }
                    }
                    for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                    {
                        //reference array
                        int rc[3]={lc[wrap%3],lc[(1+wrap)%3],lc[(2+wrap)%3]};

                        //work array
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
                                    int wc_orig[3]={wc[0],wc[1],wc[2]};
                                    //check the boundaries for each component
                                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                                    {

                                        if(wc[xyz]<0)
                                        {
                                            wc[xyz]=geom::zpdim[xyz]*geom::Nk[xyz]+wc[xyz];

                                        }


                                    }

                                    if(check(wc[0],wc[1],wc[2])==0)
                                    {
                                        if(geom::coords(wc[0],wc[1],wc[2],0)>-2)
                                        {
                                            intmat::Nrxx(wc[0],wc[1],wc[2])+=(J(i,0,0)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrxy(wc[0],wc[1],wc[2])+=(J(i,0,1)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrxz(wc[0],wc[1],wc[2])+=(J(i,0,2)/(mat::muB*mat::mustore[0]));
                                            intmat::Nryx(wc[0],wc[1],wc[2])+=(J(i,1,0)/(mat::muB*mat::mustore[0]));
                                            intmat::Nryy(wc[0],wc[1],wc[2])+=(J(i,1,1)/(mat::muB*mat::mustore[0]));
                                            intmat::Nryz(wc[0],wc[1],wc[2])+=(J(i,1,2)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrzx(wc[0],wc[1],wc[2])+=(J(i,2,0)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrzy(wc[0],wc[1],wc[2])+=(J(i,2,1)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrzz(wc[0],wc[1],wc[2])+=(J(i,2,2)/(mat::muB*mat::mustore[0]));

                                            check(wc[0],wc[1],wc[2])=1;
                                            //intmap << wc_orig[0] << "\t" << wc_orig[1] << "\t" << wc_orig[2] << "\t" << sqrt(double(wc_orig[0]*wc_orig[0]+wc_orig[1]*wc_orig[1]+wc_orig[2]*wc_orig[2]));
                                            for(unsigned int jc1 = 0 ;jc1 < 3 ;jc1++)
                                            {
                                                for(unsigned int jc2 = 0 ; jc2< 3 ; jc2++)
                                                {
                                                    //intmap << "\t" << J(i,jc1,jc2);
                                                }
                                            }
                                            //intmap << std::endl;
                                        }
                                        else
                                        {
                                            error::errPreamble(__FILE__,__LINE__);
                                            error::errMessage("That interaction does not exist");
                                        }

                                        counter++;
                                    }
                                }
                            }
                        }
                    }
                    //store the original vector
                    lc[0]=kvec(i,0);
                    lc[1]=kvec(i,1);
                    lc[2]=kvec(i,2);
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        if(geom::zpcheck==false && abs(lc[xyz])>geom::dim[xyz]*geom::Nk[xyz]/2)
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("Interactions are out of range");
                        }
                    }

                    for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                    {
                        //reference array
                        int rc[3]={lc[wrap%3],lc[(1+wrap)%3],lc[(2+wrap)%3]};
                        //work array
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
                                    int wc_orig[3]={wc[0],wc[1],wc[2]};
                                    //check the boundaries for each component
                                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                                    {
                                        if(wc[xyz]<0)
                                        {
                                            wc[xyz]=geom::zpdim[xyz]*geom::Nk[xyz]+wc[xyz];
                                        }


                                    }

                                    if(check(wc[0],wc[1],wc[2])==0)
                                    {
                                        if(geom::coords(wc[0],wc[1],wc[2],0)>-2)
                                        {
                                            intmat::Nrxx(wc[0],wc[1],wc[2])+=(J(i,0,0)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrxy(wc[0],wc[1],wc[2])+=(J(i,0,1)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrxz(wc[0],wc[1],wc[2])+=(J(i,0,2)/(mat::muB*mat::mustore[0]));
                                            intmat::Nryx(wc[0],wc[1],wc[2])+=(J(i,1,0)/(mat::muB*mat::mustore[0]));
                                            intmat::Nryy(wc[0],wc[1],wc[2])+=(J(i,1,1)/(mat::muB*mat::mustore[0]));
                                            intmat::Nryz(wc[0],wc[1],wc[2])+=(J(i,1,2)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrzx(wc[0],wc[1],wc[2])+=(J(i,2,0)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrzy(wc[0],wc[1],wc[2])+=(J(i,2,1)/(mat::muB*mat::mustore[0]));
                                            intmat::Nrzz(wc[0],wc[1],wc[2])+=(J(i,2,2)/(mat::muB*mat::mustore[0]));
                                            check(wc[0],wc[1],wc[2])=1;
                                            //intmap << wc_orig[0] << "\t" << wc_orig[1] << "\t" << wc_orig[2] << "\t" << sqrt(double(wc_orig[0]*wc_orig[0]+wc_orig[1]*wc_orig[1]+wc_orig[2]*wc_orig[2]));
                                            for(unsigned int jc1 = 0 ;jc1 < 3 ;jc1++)
                                            {
                                                for(unsigned int jc2 = 0 ; jc2< 3 ; jc2++)
                                                {
                                                    //intmap << "\t" << J(i,jc1,jc2);
                                                }
                                            }
                                            //intmap << std::endl;
                                        }
                                        else
                                        {
                                            error::errPreamble(__FILE__,__LINE__);
                                            error::errMessage("That interaction does not exist");
                                        }

                                        counter++;
                                    }
                                }
                            }
                        }
                    }
                    if(counter!=numint(i))
                    {
                        std::cerr << "Found " << counter << " should be " << numint(i) << "\tfor lc=[" << lc[0] <<" , " << lc[1] << " , " << lc[2] << "]" <<std::endl;
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("Number of interactions is not correct");
                    }
                }
                check.clear();
            }
            else if(method=="direct")
            {
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
                        if(geom::zpcheck==false && abs(oc[rc])>geom::dim[rc]*geom::Nk[rc]/2)
                        {
                            //std::cout << oc[0] << "\t" << oc[1] << "\t" << oc[2] << std::endl;
                            //std::cout << geom::dim[0]*geom::Nk[0]/2 << "\t" << geom::dim[1]*geom::Nk[1]/2 << "\t" << geom::dim[2]*geom::Nk[2]/2 << std::endl;
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("Interactions are out of range");
                        }

                        dist+=(double(c[rc])*double(c[rc]));
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
                    /*if(dist<2.0)
                      {
                      std::cout << oc[0] << "\t" << oc[1] << "\t" << oc[2] << "\t" << J[0][0] << "\t" << J[0][1] << "\t" << J[0][2] << "\t";
                      std::cout << J[1][0] << "\t" << J[1][1] << "\t" << J[1][2] << "\t";
                      std::cout << J[2][0] << "\t" << J[2][1] << "\t" << J[2][2] << std::endl;
                      }

                    //std::cout << "Jij:\n" << J[0][0] << "\t" << J[0][1] << "\t" << J[0][2] << std::endl;
                    //std::cout << J[1][0] << "\t" << J[1][1] << "\t" << J[1][2] << std::endl;
                    //std::cout << J[2][0] << "\t" << J[2][1] << "\t" << J[2][2] << std::endl;
                    //std::cin.get();*/

                    if(check(c[0],c[1],c[2])==0)//then we do not already have an interaction there
                    {

                        check(c[0],c[1],c[2])=1;
                        counter++;
                        if(geom::coords(c[0],c[1],c[2],0)>-2)
                        {
                            intmat::Nrxx(c[0],c[1],c[2])+=((J[0][0])/(mat::muB*mat::mustore[0]));
                            intmat::Nrxy(c[0],c[1],c[2])+=((J[0][1])/(mat::muB*mat::mustore[0]));
                            intmat::Nrxz(c[0],c[1],c[2])+=((J[0][2])/(mat::muB*mat::mustore[0]));
                            intmat::Nryx(c[0],c[1],c[2])+=((J[1][0])/(mat::muB*mat::mustore[0]));
                            intmat::Nryy(c[0],c[1],c[2])+=((J[1][1])/(mat::muB*mat::mustore[0]));
                            intmat::Nryz(c[0],c[1],c[2])+=((J[1][2])/(mat::muB*mat::mustore[0]));
                            intmat::Nrzx(c[0],c[1],c[2])+=((J[2][0])/(mat::muB*mat::mustore[0]));
                            intmat::Nrzy(c[0],c[1],c[2])+=((J[2][1])/(mat::muB*mat::mustore[0]));
                            intmat::Nrzz(c[0],c[1],c[2])+=((J[2][2])/(mat::muB*mat::mustore[0]));

                            //intmat::Nrzz(c[0],c[1],c[2])+=((J[2][2]+2.0*2.0*1.6e-19*1e-3)/(mat::muB*mat::mu));
                            /*                std::cout << "Interaction: " << i << "\nJij:\n" << intmat::Nrxx(c[0],c[1],c[2]) << "\t" << intmat::Nrxy(c[0],c[1],c[2]) << "\t" << intmat::Nrxz(c[0],c[1],c[2]) << std::endl;
                                              std::cout <<  intmat::Nryx(c[0],c[1],c[2]) << "\t" << intmat::Nryy(c[0],c[1],c[2]) << "\t" << intmat::Nryz(c[0],c[1],c[2]) << std::endl;
                                              std::cout <<  intmat::Nrzx(c[0],c[1],c[2]) << "\t" << intmat::Nrzy(c[0],c[1],c[2]) << "\t" << intmat::Nrzz(c[0],c[1],c[2]) << std::endl;
                                              std::cin.get();*/
                        }
                        else
                        {
                            //std::cout << oc[0] << "\t" << oc[1] << "\t" << oc[2] << "\t" <<  c[0] << "\t" << c[1] << "\t" << c[2] << std::endl;
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("You are trying to add an interaction to an empty mesh point.");
                        }
                    }
                    else
                    {
                        std::cerr << c[0] << "\t" << c[1] << "\t" << c[2] << std::endl;
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
        else
        {
            //the only reason for not using the interaction matrix should
            //be if we do not have translational invariance. Therefore
            //we don't need to zero pad
            if(geom::zpcheck==true)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You are using a neighbourlist with zero padding. Unless you have coded up multiple species then you don't need to....");
            }
            for(unsigned int pb = 0 ; pb < 3 ; pb++)
            {
                pbc[pb]=setting["pbc"][pb];
            }
            FIXOUTVEC(config::Info,"Boundary conditions:",config::isTF(pbc[0]),config::isTF(pbc[1]),config::isTF(pbc[2]));
            xadj.resize(geom::nspins+1);
            xadj.IFill(0);
            //create a temporary vector for the neighbour list
            std::vector<unsigned int> tadjncy;
            //and one for each of the isotropic exchange interactions
            std::vector<double> tJxx,tJyy,tJzz;

            //read in the exchange interactions
            if(method=="permute")
            {
                libconfig::Setting &exchset = exchcfg.lookup("exchange");
                exchset.lookupValue("Num_Shells",num_shells);
                assert(num_shells>0);
                FIXOUT(config::Info,"Reading exchange information for:" << num_shells << " shells" << std::endl);
                numint.resize(num_shells);
                exchvec.resize(num_shells,3);
                kvec.resize(num_shells,3);
                J.resize(num_shells,3,3);
                //Read the units of the exchange energy and scale appropriately to Joules
                exchset.lookupValue("units",enerType);
                FIXOUT(config::Info,"Units of exchange:" << enerType << std::endl);

                for(unsigned int i = 0 ; i < num_shells ; i++)
                {
                    std::stringstream nisstr,evsstr;
                    std::string nistr,evstr;
                    nisstr << "NumInt" << i+1;
                    nistr=nisstr.str();
                    evsstr << "Shell" << i+1 << "Vec";
                    evstr=evsstr.str();
                    exchset.lookupValue(nistr.c_str(),numint[i]);
                    config::printline(config::Info);
                    config::Info << "Shell " << i+1;
                    FIXOUT(config::Info," number of interactions:" << numint[i] << std::endl);
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {
                        exchvec(i,j)=exchset[evstr.c_str()][j];
                        kvec(i,j)=int(exchvec(i,j)*geom::Nk[j]);
                    }
                    FIXOUTVEC(config::Info,"Vectors:",exchvec(i,0),exchvec(i,1),exchvec(i,2));
                    FIXOUTVEC(config::Info,"On K-mesh:",kvec(i,0),kvec(i,1),kvec(i,2));
                    FIXOUT(config::Info,"Distance in k-points:" << sqrt(double(kvec(i,0)*kvec(i,0)+kvec(i,1)*kvec(i,1)+kvec(i,2)*kvec(i,2))));
                    config::Info << std::endl;
                    //loop over alpha (column in exchange tensor)
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {

                        std::stringstream Jsstr;
                        Jsstr << "J" << i+1 << "_" << j+1;
                        std::string Jstr;
                        Jstr=Jsstr.str();
                        //loop over beta (row in exchange tensor)
                        for(unsigned int k = 0 ; k < 3 ; k++)
                        {
                            J(i,j,k) = exchset[Jstr.c_str()][k];
                            if(enerType=="mRy")
                            {
                                J(i,j,k)*=2.179872172e-18;
                                //now to milli
                                J(i,j,k)*=1.0e-3;
                            }
                            else if(enerType=="eV")
                            {
                                J(i,j,k)*=1.602176565e-19;
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
                        FIXOUTVEC(config::Info,Jstr,J(i,j,0),J(i,j,1),J(i,j,2));
                    }

                    config::Info << std::endl;

                }
                //counter to keep track of adjncy size
                unsigned int adjncycount=0;
                FIXOUT(config::Info,"Finding neighbours...." << std::flush);
                //now we have read in the exchange matrix we need to loop over the atoms and determine the neighbours
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
//                    std::cout << "Finding neighbours for spin " << i << std::flush;
                    Array3D<int> check(geom::zpdim[0]*geom::Nk[0],geom::zpdim[1]*geom::Nk[1],geom::zpdim[2]*geom::Nk[2]);//This array is used to check if we already have the Jij for this interaction
                    check.IFill(0);

                    xadj[i]=adjncycount;
//                    std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
                    int sc[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                    //loop over the "shells" of interactions
                    for(unsigned int j = 0 ; j < num_shells ; j++)
                    {
//                    std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
                        //for debugging
                        unsigned int neighcount=0;
                        //store the original vector to permute over
                        unsigned int lc[3]={kvec(j,0),kvec(j,1),kvec(j,2)};
//                    std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
                        //permute over the interactions
                        for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                        {
//                    std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
                            //reference array
                            int rc[3]={lc[wrap%3],lc[(1+wrap)%3],lc[(2+wrap)%3]};
                            //work array
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
                                        int luc[3]={sc[0]+wc[0],sc[1]+wc[1],sc[2]+wc[2]};
//std::cout << "Before\t" << luc[0] << "\t" << luc[1] << "\t" << luc[2] << std::endl;
                                        //check boundary conditions
                                        for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                                        {
                                            if(pbc[xyz])
                                            {
                                                if(luc[xyz]<0)
                                                {
                                                    luc[xyz]=luc[xyz]%int(geom::dim[xyz])+abs(geom::dim[xyz]);
                                                }
                                                else if(luc[xyz]>=geom::dim[xyz])
                                                {
                                                    luc[xyz]=luc[xyz]%int(geom::dim[xyz]);
                                                }
                                            }
                                        }
//std::cout << "After\t" << luc[0] << "\t" << luc[1] << "\t" << luc[2] << std::endl;
                                        //if we have found a neighbour add one to the list
                                        if(luc[0] < geom::dim[0] && luc[1] < geom::dim[1] && luc[2] < geom::dim[2] && geom::coords(luc[0],luc[1],luc[2],0)>-1 && check(luc[0],luc[1],luc[2])==0)
                                        {
                        //                    std::cout << luc[0] << "\t" << luc[1] << "\t" << luc[2] << std::endl;
                                            tadjncy.push_back(1);
                                            tadjncy[adjncycount]=geom::coords(luc[0],luc[1],luc[2],0);
                                            if(geom::coords(luc[0],luc[1],luc[2],0)>geom::nspins)
                                            {
                                                error::errPreamble(__FILE__,__LINE__);
                                                error::errMessage("Neighbour larger than nspins");
                                            }
                                            tJxx.push_back(1);
                                            tJyy.push_back(1);
                                            tJzz.push_back(1);
                                            tJxx[adjncycount]=J(j,0,0)/(mat::muB*mat::mu[i]);
                                            tJyy[adjncycount]=J(j,1,1)/(mat::muB*mat::mu[i]);
                                            tJzz[adjncycount]=J(j,2,2)/(mat::muB*mat::mu[i]);
                                            adjncycount++;
                                            neighcount++;
                                            check(luc[0],luc[1],luc[2])=1;
                                        }
                                    }
                                }
                            }
                    //std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
                        }
                        if(neighcount!=numint(j) && pbc[0] && pbc[1] && pbc[2])
                        {

                            std::cerr << neighcount << "\tshould be\t" << numint(j) << std::endl;
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("Periodic boundary conditions on all size did not find correct number of neighbours");
                        }
                        neighcount=0;

                    }
//                    std::cout << "\tDone" << std::endl;
                }
                SUCCESS(config::Info);

                xadj[geom::nspins]=adjncycount;
                adjncy.resize(adjncycount);
                Jxx.resize(adjncycount);
                Jyy.resize(adjncycount);
                Jzz.resize(adjncycount);
                for(unsigned int i = 0 ; i < adjncycount ; i++)
                {
                    adjncy[i]=tadjncy[i];
                    Jxx[i]=tJxx[i];
                    Jyy[i]=tJyy[i];
                    Jzz[i]=tJzz[i];
                }
                for(unsigned int i = 0 ; i < geom::nspins; i++)
                {
                    //std::cout << "Atom " << i << " has neighbours\n";
                    for(unsigned int n=xadj[i] ; n < xadj[i+1] ; n++)
                    {
                    //    std::cout << adjncy[n] << "\t";
                    }
                    //std::cout << std::endl;
                }
            }
        }
//        exit(0);
    }
}
