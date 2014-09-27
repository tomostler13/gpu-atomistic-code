// File: exch.cpp
// Author: Tom Ostler
// Created: 18 Jan 2013
// Last-modified: 26 Sep 2014 16:36:19
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/llg.h"
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
    unsigned int max_shells=0;
    Array3D<unsigned int> numint;
    Array2D<unsigned int> shell_list;
    Array4D<double> exchvec;
    Array4D<unsigned int> kvec;
    Array5D<double> J;
    std::string enerType;
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
        numint.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells);
        exchvec.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells,3);
        kvec.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells,3);
        J.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells,3,3);
        shell_list.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS());
        numint.IFill(0);exchvec.IFill(0);kvec.IFill(0);J.IFill(0);shell_list.IFill(0);
        if(geom::nspins>1)//then really we don't want any exchange anyway.
        {
            if(method=="permute")
            {
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
                                kvec(s1,s2,i,j)=int(exchvec(s1,s2,i,j)*geom::Nk[j]);
                            }
                            FIXOUTVEC(config::Info,"Vectors:",exchvec(s1,s2,i,0),exchvec(s1,s2,i,1),exchvec(s1,s2,i,2));
                            FIXOUTVEC(config::Info,"On K-mesh:",kvec(s1,s2,i,0),kvec(s1,s2,i,1),kvec(s1,s2,i,2));
                            FIXOUT(config::Info,"Distance in k-points:" << sqrt(static_cast<double>(kvec(s1,s2,i,0)*kvec(s1,s2,i,0)+kvec(s1,s2,i,1)*kvec(s1,s2,i,1)+kvec(s1,s2,i,2)*kvec(s1,s2,i,2))));
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
                            lc[0]=kvec(s1,s2,i,0);
                            lc[1]=kvec(s1,s2,i,1);
                            lc[2]=kvec(s1,s2,i,2);
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
                                            //This array is purely for outputting the exchange information to the log file
                                            int owc[3]={wc[0],wc[1],wc[2]};
                                            bool checkmonolayer[3]={false,false,false};
                                            //check the boundaries for each component
                                            for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                                            {
                                                if(wc[xyz]<0)
                                                {
                                                    if(geom::dim[xyz]*geom::Nk[xyz] < 2 && wc[xyz]<0)
                                                    {
                                                        checkmonolayer[xyz]=true;
                                                    }
                                                    wc[xyz]=geom::zpdim[xyz]*geom::Nk[xyz]+wc[xyz];
                                                }
                                            }//end of xyz loop
                                            if(check(wc[0],wc[1],wc[2])==0)
                                            {
                                                config::Log << "Interaction Vector:  [" << owc[0] << "," << owc[1] << "," << owc[2] << "]\t -> [" << wc[0] << "," << wc[1] << "," << wc[2] << "]" << std::endl;
                                                //loop over the elements of the interaction tensor
                                                for(unsigned int alpha = 0 ; alpha < 3 ; alpha++)
                                                {
                                                    for(unsigned int beta = 0 ; beta < 3 ; beta++)
                                                    {
                                                        if(checkmonolayer[0]==true || checkmonolayer[1]==true || checkmonolayer[2]==true)
                                                        {
                                                            //do nothing
                                                        }
                                                        else
                                                        {
                                                            intmat::Nrab(s1,s2,alpha,beta,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,alpha,beta)/(geom::ucm.GetMuBase(s1)*llg::muB));
                                                        }
                                                    }//end of beta loop
                                                }//end of alpha loop

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
