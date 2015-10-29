// File: exch_interaction_matrix.cpp
// Author: Tom Ostler
// Created: 05 Dec 2014
// Last-modified: 29 Oct 2015 10:02:51
// This source file was added to tidy up the file exch.cpp
// because it was becoming cumbersome to work with. This
// source file calculates the interaction matrices based
// on the exchange values that are read in (see
// exch_determine_exchange_matrix.cpp
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
    void fft()
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
                    config::Log << "Shell " << i << " of " << shell_list(s1,s2) << " -> " << numint(s1,s2,i) << " interactions" << std::endl;
                    unsigned int counter=0;
                    int lc[3]={0,0,0};
                    lc[0]=static_cast<unsigned int>(exchvec(s1,s2,i,0)*static_cast<double>(geom::Nk[0])*evs[0]+0.5);
                    lc[1]=static_cast<unsigned int>(exchvec(s1,s2,i,1)*static_cast<double>(geom::Nk[1])*evs[1]+0.5);
                    lc[2]=static_cast<unsigned int>(exchvec(s1,s2,i,2)*static_cast<double>(geom::Nk[2])*evs[2]+0.5);
                    for(unsigned int wrap = 0 ; wrap < 3 ; wrap++)
                    {
                        //reference array
                        int rc[3]={lc[wrap%3],lc[(1+wrap)%3],lc[(2+wrap)%3]};
                        //swap the last two elements to get all permutations
                        for(unsigned int swap = 0 ; swap < 2 ; swap++)
                        {
                            if(swap==1)
                            {
                                int el2=rc[1];
                                int el3=rc[2];
                                rc[2]=el2;
                                rc[1]=el3;
                            }

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
                                                std::cerr << __FILE__ << "\t" << __LINE__ << "\tThis section of code is wrong." << std::endl;
                                                exit(0);
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
                        }//end of swap loop
                    }//end of wrap loop
                    if(counter!=numint(s1,s2,i))
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream sstr;
                        sstr << "The number of interactions between species " << s1 << " and " << s2 << " in shell " << i << "(of " << shell_list(s1,s2) << ") should be " << numint(s1,s2,i) << ", instead the number counted was " << counter;
                        std::string str=sstr.str();
                        error::errMessage(str.c_str());
                    }//end of counter check if statement
                    counter=0;
                    config::printline(config::Log);
                }//end of shell list loop
            }//end of s2 loop
        }//end of s1 loop
    }
    void directfft()
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
//                    config::Log << "Shell " << shell_list(s1,s2) << " of " << numint(s1,s2,i) << std::endl;
                    unsigned int counter=0;
                    int lc[3]={0,0,0};
                    //std::cout << s1 << "\t" << s2 << "\t" << i << "\t" << exchvec(s1,s2,i,0) << std::endl;
                    //std::cout << evs[0] << std::endl;
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
                                            intmat::Nrab(s1,s2,0,1,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,1,0))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                            // Nxz = 1/2(Jxz-Jzx)
                                            intmat::Nrab(s1,s2,0,2,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,0,2))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                            // Nyx = 1/2(Jxy-Jyx)
                                            intmat::Nrab(s1,s2,1,0,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,1,0))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                            // Nyz = 1/2(Jzy-Jyz)
                                            intmat::Nrab(s1,s2,1,2,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,1,2))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                            // Nzx = 1/2(Jzx - Jxz)
                                            intmat::Nrab(s1,s2,2,0,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,2,0))/(geom::ucm.GetMuBase(s1)*llg::muB);
                                            // Nzy = 1/2(Jyz-Jzy)
                                            intmat::Nrab(s1,s2,2,1,wc[0],wc[1],wc[2])[0]+=(J(s1,s2,i,2,1))/(geom::ucm.GetMuBase(s1)*llg::muB);

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
/*                    if(counter!=numint(s1,s2,i))
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream sstr;
                        sstr << "The number of interactions between species " << s1 << " and " << s2 << " in shell " << shell_list(s1,s2) << " should be " << numint(s1,s2,i) << ", instead the number counted was " << counter;
                        std::string str=sstr.str();
                        error::errMessage(str.c_str());
                    }//end of counter check if statement
                    */
                    counter=0;
                    config::printline(config::Log);
                }//end of shell list loop
            }//end of s2 loop
        }//end of s1 loop
    }
}
