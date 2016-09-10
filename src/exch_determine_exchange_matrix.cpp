// File: exch_determine_exchange_matrix.cpp
// Author: Tom Ostler
// Created: 05 Dec 2014
// Last-modified: 10 Sep 2016 13:38:32
// This source file was added to tidy up the file exch.cpp
// because it was becoming cumbersome to work with. This
// source file calculates the CSR neighbourlist
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
    unsigned int maxNoSpecInt;
    void get_exch_unitcell(int argc,char *argv[])
    {
        FIXOUT(config::Info,"Determining interactions within and between unit cells" << std::endl);
        J.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_int,3,3);
        for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
        {
            for(unsigned int s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
            {
                config::printline(config::Info);
                FIXOUT(config::Info,"Exchange interaction between species:" << s1 << " and " << s2 << std::endl);
                std::stringstream sstr_int;
                sstr_int << "exchange" << "_" << s1 << "_" << s2;
                std::string str_int = sstr_int.str();
                libconfig::Setting &exchset = exchcfg.lookup(str_int.c_str());

                if(exchset.lookupValue("Num_Interactions",shell_list(s1,s2)))
                {
                    FIXOUT(config::Info,"Reading information for:" << shell_list(s1,s2) << " interactions" << std::endl);
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Could not read the number of interactions between species " << s1 << " and " << s2;
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
                if(exchset.lookupValue("units",enerType))
                {
                    FIXOUT(config::Info,"Units:" << enerType << std::endl);
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Error reading energy type, check your exchange file.");
                }

                for(unsigned int i = 0 ; i < shell_list(s1,s2) ; i++)
                {
                    std::stringstream evsstr;
                    std::string evstr;

                    evsstr << "UnitCell" << i+1;
                    evstr=evsstr.str();
                    //get the unit cell for interaction
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {
                        try
                        {
                            exchvec(s1,s2,i,j)=exchset[evstr.c_str()][j];
                        }
                        catch(const libconfig::SettingNotFoundException &snf)
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            std::stringstream errsstr;
                            errsstr << "Could not read vector for interaction number, " << i+1 << " check your exchange file.";
                            std::string errstr=errsstr.str();
                            error::errMessage(errstr);
                        }
                    }
                    FIXOUTVEC(config::Info,"Unit cell vector:",exchvec(s1,s2,i,0),exchvec(s1,s2,i,1),exchvec(s1,s2,i,2));
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {
                        std::stringstream Jsstr;
                        Jsstr << "J" << i+1 << "_" << j+1;
                        std::string Jstr;
                        Jstr=Jsstr.str();
                        for(unsigned int k = 0 ; k < 3 ; k++)
                        {
                            try
                            {
                                J(s1,s2,i,j,k) = exchset[Jstr.c_str()][k];
                            }
                            catch(const libconfig::SettingNotFoundException &snf)
                            {
                                std::stringstream errsstr;
                                errsstr << "Could not exchange interaction number " << i+1 << " tensor components " << j << " and " << k << ", for interaction between species " << s1 << " and " << s2 << " (" << snf.getPath() << ")";
                                std::string errstr=errsstr.str();
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage(errstr);
                            }
                            if(enerType=="mRy")
                            {
                                J(s1,s2,i,j,k)*=2.179872172e-18; //now to milli
                                J(s1,s2,i,j,k)*=1.0e-3;
                            }
                            else if(enerType=="eV")
                            {
                                J(s1,s2,i,j,k)*=1.602176565e-19;
                            }
                            else if(enerType=="meV")
                            {
                                J(s1,s2,i,j,k)*=1.602176565e-22;
                            }
                            else if(enerType=="J" || enerType=="Joules" || enerType=="joules")
                            {
                                //do nothing
                            }
                            else if(enerType=="Ry")
                            {
                                J(s1,s2,i,j,k)*=2.179872172e-18;
                            }
                            else
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("Units of exchange energy not recognised");
                            }

                        }
                        FIXOUTVEC(config::Info,Jstr,J(s1,s2,i,j,0),J(s1,s2,i,j,1),J(s1,s2,i,j,2));
                    }
                }
            }
        }
    }
    void get_exch_permute(int argc,char *argv[])
    {
        //first read the exchange constants
        J.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells,3,3);
        Array2D<unsigned int> maxNoInt;
        maxNoInt.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS());maxNoInt.IFill(0);
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
                    maxNoInt(s1,s2)+=numint(s1,s2,i);
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
                            else if(enerType=="Ry")
                            {
                                 J(s1,s2,i,j,k)*=2.179872172e-18;
                            }
                            else if(enerType=="eV")
                            {
                                J(s1,s2,i,j,k)*=1.602176565e-19;
                            }
                            else if(enerType=="meV")
                            {
                                J(s1,s2,i,j,k)*=1.602176565e-19*1e-3;
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
                if(maxNoInt(s1,s2)>maxNoSpecInt)
                {
                    maxNoSpecInt=maxNoInt(s1,s2);
                }
            }
        }
        FIXOUT(config::Info,"Max number of interactions a given spin can have is:" << maxNoSpecInt << std::endl);
    }

    void get_exch_mapint(int argc,char *argv[])
    {
        FIXOUT(config::Info,"Determining interactions from a list of integer mesh lookups" << std::endl);
        //first read the exchange constants
        J.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_int,3,3);
        for(unsigned int s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
        {
            for(unsigned int s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
            {
                config::printline(config::Info);
                FIXOUT(config::Info,"Exchange interaction between species:" << s1 << " and " << s2 << std::endl);
                std::stringstream sstr_int;
                sstr_int << "exchange" << "_" << s1 << "_" << s2;
                std::string str_int = sstr_int.str();
                libconfig::Setting &exchset = exchcfg.lookup(str_int.c_str());

                //If the interactions are "direct" then the array shell_list actually stores the total number
                //of interactions for the interaction between species s1 and s2
                exchset.lookupValue("Num_Interactions",shell_list(s1,s2));

                FIXOUT(config::Info,"Reading information for:" << shell_list(s1,s2) << " interactions" << std::endl);

                exchset.lookupValue("units",enerType);
                for(unsigned int i = 0 ; i < shell_list(s1,s2) ; i++)
                {
                    std::stringstream evsstr;
                    std::string evstr;

                    evsstr << "Interaction" << i+1 << "Vec";
                    evstr=evsstr.str();
                    //get the vector for interaction i
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {
                        exchvec(s1,s2,i,j)=exchset[evstr.c_str()][j];

                    }
                    config::Info << "Interaction " << i << " ";
                    FIXOUTVEC(config::Info," vector:",exchvec(s1,s2,i,0),exchvec(s1,s2,i,1),exchvec(s1,s2,i,2));
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
                            else if(enerType=="meV")
                            {
                                J(s1,s2,i,j,k)*=1.602176565e-22;
                            }
                            else if(enerType=="J" || enerType=="Joules" || enerType=="joules")
                            {
                                //do nothing
                            }
                            else if(enerType=="Ry")
                            {
                                J(s1,s2,i,j,k)*=2.179872172e-18;
                            }
                            else
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("Units of exchange energy not recognised");
                            }

                        }
                        FIXOUTVEC(config::Info,Jstr,J(s1,s2,i,j,0),J(s1,s2,i,j,1),J(s1,s2,i,j,2));
                    }
                }
            }
        }
    }
    void hybrid()
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

                lc[0]=static_cast<unsigned int>(exchvec(s1,s1,i,0)*static_cast<double>(geom::Nk[0])*evs[0]+0.5);
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
                                        intmat::hNrab(0,1,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nxz = 1/2(Jxz-Jzx)
                                        intmat::hNrab(0,2,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nyx = 1/2(Jxy-Jyx)
                                        intmat::hNrab(1,0,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nyz = 1/2(Jzy-Jyz)
                                        intmat::hNrab(1,2,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nzx = 1/2(Jzx - Jxz)
                                        intmat::hNrab(2,0,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nzy = 1/2(Jyz-Jzy)
                                        intmat::hNrab(2,1,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);

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



    }
    void directhybrid()
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
                                        intmat::hNrab(0,1,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nxz = 1/2(Jxz-Jzx)
                                        intmat::hNrab(0,2,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nyx = 1/2(Jxy-Jyx)
                                        intmat::hNrab(1,0,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nyz = 1/2(Jzy-Jyz)
                                        intmat::hNrab(1,2,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nzx = 1/2(Jzx - Jxz)
                                        intmat::hNrab(2,0,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);
                                        // Nzy = 1/2(Jyz-Jzy)
                                        intmat::hNrab(2,1,plane,wc[1],wc[2])[0]/=(geom::ucm.GetMuBase(s1)*llg::muB);

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



    }
}
