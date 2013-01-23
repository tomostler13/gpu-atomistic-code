// File: exch.cpp
// Author: Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 23 Jan 2013 11:52:26
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/mat.h"
#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

namespace exch
{
    unsigned int num_shells = 0;
    Array<unsigned int> numint;
    Array2D<double> exchvec;
    Array2D<unsigned int> kvec;
    Array3D<double> J;
    std::string enerType;
    void initExch(int argc,char *argv[])
    {
        std::string readMethod,readFile;
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
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Method of reading exchange file not recognised.");
        }
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
                                    intmat::Nrxx(wc[0],wc[1],wc[2])+=(J(i,0,0)/(mat::muB*mat::mu));
                                    intmat::Nrxy(wc[0],wc[1],wc[2])+=(J(i,0,1)/(mat::muB*mat::mu));
                                    intmat::Nrxz(wc[0],wc[1],wc[2])+=(J(i,0,2)/(mat::muB*mat::mu));
                                    intmat::Nryx(wc[0],wc[1],wc[2])+=(J(i,1,0)/(mat::muB*mat::mu));
                                    intmat::Nryy(wc[0],wc[1],wc[2])+=(J(i,1,1)/(mat::muB*mat::mu));
                                    intmat::Nryz(wc[0],wc[1],wc[2])+=(J(i,1,2)/(mat::muB*mat::mu));
                                    intmat::Nrzx(wc[0],wc[1],wc[2])+=(J(i,2,0)/(mat::muB*mat::mu));
                                    intmat::Nrzy(wc[0],wc[1],wc[2])+=(J(i,2,1)/(mat::muB*mat::mu));
                                    intmat::Nrzz(wc[0],wc[1],wc[2])+=(J(i,2,2)/(mat::muB*mat::mu));
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
            lc[0]=kvec(i,2);
            lc[1]=kvec(i,1);
            lc[2]=kvec(i,0);
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
                                    intmat::Nrxx(wc[0],wc[1],wc[2])+=(J(i,0,0)/(mat::muB*mat::mu));
                                    intmat::Nrxy(wc[0],wc[1],wc[2])+=(J(i,0,1)/(mat::muB*mat::mu));
                                    intmat::Nrxz(wc[0],wc[1],wc[2])+=(J(i,0,2)/(mat::muB*mat::mu));
                                    intmat::Nryx(wc[0],wc[1],wc[2])+=(J(i,1,0)/(mat::muB*mat::mu));
                                    intmat::Nryy(wc[0],wc[1],wc[2])+=(J(i,1,1)/(mat::muB*mat::mu));
                                    intmat::Nryz(wc[0],wc[1],wc[2])+=(J(i,1,2)/(mat::muB*mat::mu));
                                    intmat::Nrzx(wc[0],wc[1],wc[2])+=(J(i,2,0)/(mat::muB*mat::mu));
                                    intmat::Nrzy(wc[0],wc[1],wc[2])+=(J(i,2,1)/(mat::muB*mat::mu));
                                    intmat::Nrzz(wc[0],wc[1],wc[2])+=(J(i,2,2)/(mat::muB*mat::mu));
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
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Number of interactions is not correct");
            }
        }
        check.clear();
    }
}
