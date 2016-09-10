// File: util_mag.cpp
// Author:Tom Ostler
// Created: 15 Dec 2014
// Last-modified: 10 Sep 2016 18:56:57
// Contains useful functions and classes
// that pertain to magnetization
#include "../inc/util.h"
#include "../inc/llg.h"
#include "../inc/arrays.h"
#include "../inc/fields.h"
#include "../inc/defines.h"
#include "../inc/spins.h"
#include "../inc/intmat.h"
#include "../inc/geom.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include <string>
#include <sstream>
namespace util
{
    void calc_mag()
    {
        //This section of code can be used to add new ways to calculate the magnetization.
        //For example if one wanted to calculate the height resolved magnetization then it should
        //be added here. It should be added to the list below just incase I ever want to write a
        //manual (insert ridiculous laugh here). It may be possible that you have to create a
        //new set of arrays for storing and calculating the magnetization. Using the above example
        //(m(height)) you would have to create an array in src/spins.cpp and also declare an
        //extern array in inc/spins.h.
        //
        //
        // List of arguements and what they do
        //
        // 0 - calculate the sublattice resolved magnetization
        // 1 - output the magnetization as a function of x (mu_b)
        // 2 - output the magnetization as a function of y (mu_b)
        // 3 - output the magnetization as a function of z (mu_b)
        // 4 - output the entire spin map at the spins::update frequency
        // 5 - output the layer resolved magnetization as a function of x (reduced)
        // 6 - output the layer resolved magnetization as a function of y (reduced)
        // 7 - output the layer resolved magnetization as a function of z (reduced)
        // 8 - Neel vector for a checkerboard AFM
        // 9 - Average over a number of cells in space (x,y,z)
        if(spins::mag_calc_method==0 || spins::output_mag==true)
        {
            spins::mag.IFill(0);
            //add up the magnetization for each sublattice
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                unsigned int sl = geom::sublattice[i];
                spins::mag(sl,0)+=spins::Sx[i];
                spins::mag(sl,1)+=spins::Sy[i];
                spins::mag(sl,2)+=spins::Sz[i];
            }
            //divide by the number of spins in each sublattice
            for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
            {
                double oones=1./static_cast<double>(geom::ucm.GetNES(i));
                spins::mag(i,0)*=oones;
                spins::mag(i,1)*=oones;
                spins::mag(i,2)*=oones;
            }
        }
        if(spins::mag_calc_method>0)
        {
            if(spins::mag_calc_method==1)
            {
                magx.IFill(0);
                //calculate the total magnetization in the layer
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    //check which plane the spin belongs to
                    unsigned int klu=geom::lu(i,0);
                    magx(klu,0)+=spins::Sx[i]*geom::mu[i];
                    magx(klu,1)+=spins::Sy[i]*geom::mu[i];
                    magx(klu,2)+=spins::Sz[i]*geom::mu[i];
                }

            }
            else if(spins::mag_calc_method==2)
            {
                magy.IFill(0);
                //calculate the total magnetization in the layer
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    //check which plane the spin belongs to
                    unsigned int klu=geom::lu(i,1);
                    magy(klu,0)+=spins::Sx[i]*geom::mu[i];
                    magy(klu,1)+=spins::Sy[i]*geom::mu[i];
                    magy(klu,2)+=spins::Sz[i]*geom::mu[i];
                }

            }
            else if(spins::mag_calc_method==3)
            {
                magz.IFill(0);
                //calculate the total magnetization in the layer
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    //check which plane the spin belongs to
                    unsigned int klu=geom::lu(i,2);
                    magz(klu,0)+=spins::Sx[i]*geom::mu[i];
                    magz(klu,1)+=spins::Sy[i]*geom::mu[i];
                    magz(klu,2)+=spins::Sz[i]*geom::mu[i];
                }
            }
            else if(spins::mag_calc_method==4)
            {
                // do not need to do anything, everything is done
                // when output_mag is called.
            }
            else if(spins::mag_calc_method==5)
            {
                mag_species_x.IFill(0);
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    //plane lookup
                    unsigned int klu=geom::lu(i,0);
                    //species lookup
                    unsigned int splu=geom::lu(i,3);
                    mag_species_x(klu,splu,0)+=spins::Sx[i];
                    mag_species_x(klu,splu,1)+=spins::Sy[i];
                    mag_species_x(klu,splu,2)+=spins::Sz[i];
                }
            }
            else if(spins::mag_calc_method==6)
            {
                mag_species_y.IFill(0);
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    //plane lookup
                    unsigned int klu=geom::lu(i,1);
                    //species lookup
                    unsigned int splu=geom::lu(i,3);
                    mag_species_y(klu,splu,0)+=spins::Sx[i];
                    mag_species_y(klu,splu,1)+=spins::Sy[i];
                    mag_species_y(klu,splu,2)+=spins::Sz[i];
                }
            }
            else if(spins::mag_calc_method==7)
            {
                mag_species_z.IFill(0);
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    //plane lookup
                    unsigned int klu=geom::lu(i,2);
                    //species lookup
                    unsigned int splu=geom::lu(i,3);
                    mag_species_z(klu,splu,0)+=spins::Sx[i];
                    mag_species_z(klu,splu,1)+=spins::Sy[i];
                    mag_species_z(klu,splu,2)+=spins::Sz[i];
                }
            }
            else if(spins::mag_calc_method==8)
            {
                lx2=0;ly2=0;lz2=0;
                lx1=0;ly1=0;lz2=0;
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    int mycoords[3]={geom::lu(i,0),geom::lu(i,1),geom::lu(i,2)};
                    if((mycoords[0]+mycoords[1]+mycoords[2])%2==0)
                    {//then we are on sublattice A
                        lx1+=spins::Sx[i];
                        ly1+=spins::Sy[i];
                        lz1+=spins::Sz[i];
                    }
                    else
                    {
                        lx2+=spins::Sx[i];
                        ly2+=spins::Sy[i];
                        lz2+=spins::Sz[i];
                    }
                }
                lx1=lx1/(static_cast<double>(geom::nspins)*0.5);
                ly1=ly1/(static_cast<double>(geom::nspins)*0.5);
                lz1=lz1/(static_cast<double>(geom::nspins)*0.5);
                lx2=lx2/(static_cast<double>(geom::nspins)*0.5);
                ly2=ly2/(static_cast<double>(geom::nspins)*0.5);
                lz2=lz2/(static_cast<double>(geom::nspins)*0.5);
            }
            else if(spins::mag_calc_method==9)
            {
                magdisc.IFill(0);
                //double summz=0.0;
                //double sumdmz=0.0;
                for(unsigned int i = 0 ; i < geom::nspins ; i++)
                {
                    int kx=geom::lu(i,0);
                    int ky=geom::lu(i,1);
                    int kz=geom::lu(i,2);
                    int cx=static_cast<int>(static_cast<double>(kx)/static_cast<double>(magDiscSize[0]));
                    int cy=static_cast<int>(static_cast<double>(ky)/static_cast<double>(magDiscSize[1]));
                    int cz=static_cast<int>(static_cast<double>(kz)/static_cast<double>(magDiscSize[2]));
                    magdisc(cx,cy,cz,0)+=spins::Sx[i];
                    magdisc(cx,cy,cz,1)+=spins::Sy[i];
                    magdisc(cx,cy,cz,2)+=spins::Sz[i];
                //    summz+=spins::Sz[i];
                }
                //std::cout << summz << std::endl;
                for(unsigned int cx = 0 ; cx < maxcx+1 ; cx++)
                {
                    for(unsigned int cy = 0 ; cy < maxcy+1 ; cy++)
                    {
                        for(unsigned int cz = 0 ; cz < maxcz+1 ; cz++)
                        {
                            const double oonspc=1./(static_cast<double>(nd(cx,cy,cz)));
                            for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                            {
                                magdisc(cx,cy,cz,xyz)*=oonspc;
                            }

                            //sumdmz+=magdisc(cx,cy,cz,2);
                        }
                    }
                }
                //std::cout << sumdmz << std::endl;
                //std::cin.get();


            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("The method for calculating the magnetization is not recognised.");
            }
        }
    }
    void init_output()
    {
        if(spins::mag_calc_method==0 || spins::output_mag==true)
        {
            ofs.open("mag.dat");
            if(!ofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file.");
            }
            ofs << "#File description: this file contains the magnetization as a function of time throughout the simulation." << std::endl;
            ofs << "#time - mx_0 - my_0 - mz_0 - mx_1 ...." << std::endl;
        }
        if(spins::mag_calc_method==1)
        {

            sofs.open("mag_x.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file mag_x.dat.");
            }
            sofs << "#File description: this file contains the magnetization as a function of time and x-direction in mu_B. The data is arranged in indices, one for each time." << std::endl;
            sofs << "#time - x - mx - my - mz ...." << std::endl;
            magx.resize(geom::Nk[0]*geom::dim[0],3);
        }
        else if(spins::mag_calc_method==2)
        {
            sofs.open("mag_y.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file mag_y.dat.");
            }
            sofs << "#File description: this file contains the magnetization as a function of time and y-direction in mu_B. The data is arranged in indices, one for each time." << std::endl;
            sofs << "#time - y - mx - my - mz ...." << std::endl;
            magy.resize(geom::Nk[1]*geom::dim[1],3);
        }
        else if(spins::mag_calc_method==3)
        {
            sofs.open("mag_z.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file mag_z.dat.");
            }
            sofs << "#File description: this file contains the magnetization as a function of time and z-direction in mu_B. The data is arranged in indices, one for each time." << std::endl;
            sofs << "#time - z - mx - my - mz ...." << std::endl;
            magz.resize(geom::Nk[2]*geom::dim[2],3);
        }
        else if(spins::mag_calc_method==4)
        {
            sofs.open("spec_map.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open species map file (spec_map.dat)");
            }
            sofs << "#File Description: contains the species number for each occupied lattice point. This is so that the positions do not have to be output each time for spinmap output." << std::endl;
            sofs << "# x [A] - y [A] - z [A] - atom numner - species number" << std::endl;
            for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                    {
                        if(geom::coords(i,j,k,0)>-1)
                        {
                            sofs << static_cast<double>(i)*geom::abc[0]/static_cast<double>(geom::Nk[0]) << "\t";
                            sofs << static_cast<double>(j)*geom::abc[1]/static_cast<double>(geom::Nk[1]) << "\t";
                            sofs << static_cast<double>(k)*geom::abc[2]/static_cast<double>(geom::Nk[2]) << "\t";
                            sofs << geom::coords(i,j,k,0) << "\t" << geom::coords(i,j,k,1) << std::endl;
                        }
                    }
                }
            }


            sofs.close();
        }
        else if(spins::mag_calc_method==5)
        {
            sofs.open("mag_species_x.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file mag_species_x.dat.");
            }
            sofs << "#File description: this file contains the magnetization for each species as a function of time and x-direction in mu_B. The data is arranged in indices, one for each time." << std::endl;
            sofs << "#time - x - mx(0) - my(0) - mz(0) ...." << std::endl;
            mag_species_x.resize(geom::Nk[0]*geom::dim[0],geom::ucm.GetNMS(),3);
            nspl.resize(geom::Nk[0]*geom::dim[0],geom::ucm.GetNMS());
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                //plane lookup
                unsigned int klu=geom::lu(i,0);
                //species lookup
                unsigned int splu=geom::lu(i,3);
                nspl(klu,splu)++;
            }
            for(unsigned int i = 0 ; i < geom::Nk[0]*geom::dim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::ucm.GetNMS() ; j++)
                {
                    if(nspl(i,j)<1e-50)
                    {
                        nspl(i,j)=0.0;
                    }
                    else
                    {
                        nspl(i,j)=(1./nspl(i,j));
                    }
                }
            }
        }
        else if(spins::mag_calc_method==6)
        {
            sofs.open("mag_species_y.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file mag_species_y.dat.");
            }
            sofs << "#File description: this file contains the magnetization for each species as a function of time and y-direction in mu_B. The data is arranged in indices, one for each time." << std::endl;
            sofs << "#time - y - mx(0) - my(0) - mz(0) ...." << std::endl;
            mag_species_y.resize(geom::Nk[1]*geom::dim[1],geom::ucm.GetNMS(),3);
            nspl.resize(geom::Nk[1]*geom::dim[1],geom::ucm.GetNMS());
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                //plane lookup
                unsigned int klu=geom::lu(i,1);
                //species lookup
                unsigned int splu=geom::lu(i,3);
                nspl(klu,splu)++;
            }
            for(unsigned int i = 0 ; i < geom::Nk[1]*geom::dim[1] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::ucm.GetNMS() ; j++)
                {
                    if(nspl(i,j)<1e-50)
                    {
                        nspl(i,j)=0.0;
                    }
                    else
                    {
                        nspl(i,j)=(1./nspl(i,j));
                    }
                }
            }
        }
        else if(spins::mag_calc_method==7)
        {
            sofs.open("mag_species_z.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file mag_species_z.dat.");
            }
            sofs << "#File description: this file contains the magnetization for each species as a function of time and z-direction in mu_B. The data is arranged in indices, one for each time." << std::endl;
            sofs << "#time - z - mx(0) - my(0) - mz(0) ...." << std::endl;
            mag_species_z.resize(geom::Nk[2]*geom::dim[2],geom::ucm.GetNMS(),3);
            nspl.resize(geom::Nk[2]*geom::dim[2],geom::ucm.GetNMS());
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                //plane lookup
                unsigned int klu=geom::lu(i,2);
                //species lookup
                unsigned int splu=geom::lu(i,3);
                nspl(klu,splu)++;
            }
            for(unsigned int i = 0 ; i < geom::Nk[2]*geom::dim[2] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::ucm.GetNMS() ; j++)
                {
                    if(nspl(i,j)<1e-50)
                    {
                        nspl(i,j)=0.0;
                    }
                    else
                    {
                        nspl(i,j)=(1./nspl(i,j));
                    }
                }
            }
        }
        else if(spins::mag_calc_method==8)
        {
            sofs.open("afm_check_sublat.dat");
            if(!sofs.is_open())
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not open magnetization file afm_check_sublat.dat.");
            }
            sofs << "#File description: this file contains the magnetization for each sublattice of a 3D checker board antiferromagnet." << std::endl;
            sofs << "#time - z - mx_A - my_A - mz_A - mx_B - my_B - mz_B" << std::endl;

        }
        else if(spins::mag_calc_method==9)
        {
            if(disOutForm=="single")
            {
                sofs.open("disc_mag.dat");
                if(!sofs.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open magnetization file for discretising magnetization (disc_mag.dat)\nMostly likely the disk is full or write permissions are not possible.");
                }
            }
            maxcx=0,maxcy=0,maxcz=0;

            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                int kx=geom::lu(i,0);
                int ky=geom::lu(i,1);
                int kz=geom::lu(i,2);
                int cx=static_cast<int>(static_cast<double>(kx)/static_cast<double>(magDiscSize[0]));
                int cy=static_cast<int>(static_cast<double>(ky)/static_cast<double>(magDiscSize[1]));
                int cz=static_cast<int>(static_cast<double>(kz)/static_cast<double>(magDiscSize[2]));
                if(cx>static_cast<int>(maxcx)){maxcx=cx;}
                if(cy>static_cast<int>(maxcy)){maxcy=cy;}
                if(cz>static_cast<int>(maxcz)){maxcz=cz;}
            }
            magdisc.resize(maxcx+1,maxcy+1,maxcz+1,3);
            magdisc.IFill(0);
            nd.resize(maxcx+1,maxcy+1,maxcz+1);
            nd.IFill(0);
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                int kx=geom::lu(i,0);
                int ky=geom::lu(i,1);
                int kz=geom::lu(i,2);
                int cx=static_cast<int>(static_cast<double>(kx)/static_cast<double>(magDiscSize[0]));
                int cy=static_cast<int>(static_cast<double>(ky)/static_cast<double>(magDiscSize[1]));
                int cz=static_cast<int>(static_cast<double>(kz)/static_cast<double>(magDiscSize[2]));
                nd(cx,cy,cz)++;
            }
            FIXOUTVEC(config::Info,"Number of cells in each direction for disc mag:",maxcx+1,maxcy+1,maxcz+1);
/*            for(unsigned int i = 0 ; i < maxcx+1 ; i++)
            {
                for(unsigned int j = 0 ; j < maxcy+1 ; j++)
                {
                    for(unsigned int k = 0 ; k < maxcz+1 ; k++)
                    {
                        std::cout << i << "\t" << j << "\t" << k << "\t" << nd(i,j,k) << std::endl;
                    }
                }
            }*/
        }
    }
    void output_mag(unsigned int t)
    {
        //This section of code outputs the code in a method consistent
        //with the method specified by spins::mag_calc_method.
        //If you want to introduce a new method of outputting the method
        //that is not the same as calculating the magnetization or if you
        //want to output in two ways then a new control variable should be
        //defined and read in somewhere.
        //
        //
        // List of arguements and what they do
        //
        // 0 - output the sublattice resolved magnetization
        // 1 - output the magnetization as a function of x (mu_b)
        // 2 - output the magnetization as a function of y (mu_b)
        // 3 - output the magnetization as a function of z (mu_b)
        // 4 - output the entire spin map at the spins::update frequency
        if(spins::mag_calc_method==0 || spins::output_mag==true)
        {
            ofs << static_cast<double>(t)*llg::dt << "\t";
            for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
            {
                ofs << spins::mag(s,0) << "\t" << spins::mag(s,1) << "\t" << spins::mag(s,2) << "\t";
            }
            ofs << std::endl;
        }
        if(spins::mag_calc_method==1)//along x
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[0]*geom::dim[0] ; i++)
            {
                sofs << timeid << "\t" << magx(i,0) << "\t" << magx(i,1) << "\t" << magx(i,2) << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==2)//along y
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[1]*geom::dim[1] ; i++)
            {
                sofs << timeid << "\t" << magy(i,0) << "\t" << magy(i,1) << "\t" << magy(i,2) << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==3)//along z
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[2]*geom::dim[2] ; i++)
            {
                sofs << timeid << "\t" << i << "\t" << magz(i,0) << "\t" << magz(i,1) << "\t" << magz(i,2) << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==4 && spins::mapout)
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            std::stringstream sstr;
            sstr << "spinmap_" << timeid << ".dat";
            std::string str=sstr.str();
            std::ofstream smofs(str.c_str());
            if(!smofs.is_open())
            {
                std::stringstream sstr1;
                sstr1 << "Could not open file spinmap_" << timeid << ".dat";
                std::string str1=sstr1.str();
                error::errPreamble(__FILE__,__LINE__);
                error::errWarning(str1);
            }
            for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                    {
                        int sn=geom::coords(i,j,k,0);
                        if(geom::coords(i,j,k,0)>-1)
                        {
                            smofs << spins::Sx[sn] << "\t";
                            smofs << spins::Sy[sn] << "\t";
                            smofs << spins::Sz[sn] << std::endl;
                        }
                    }
                }
            }
            smofs.close();
            if(smofs.is_open())
            {
                std::stringstream sstr1;
                sstr1 << "Could not close file spinmap_" << timeid << ".dat";
                std::string str1=sstr1.str();
                error::errPreamble(__FILE__,__LINE__);
                error::errWarning(str1);
            }

        }
        else if(spins::mag_calc_method==1)//along x
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[0]*geom::dim[0] ; i++)
            {
                const double timeid=static_cast<double>(t)*llg::dt;
                //output in index format for plotting with gnuplot
                for(unsigned int i = 0 ; i < geom::Nk[0]*geom::dim[0] ; i++)
                {
                    sofs << timeid << "\t" << magx(i,0) << "\t" << magx(i,1) << "\t" << magx(i,2) << std::endl;
                }
                sofs << std::endl << std::endl;
            }
        }
        else if(spins::mag_calc_method==2)//along y
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[1]*geom::dim[1] ; i++)
            {
                sofs << timeid << "\t" << magy(i,0) << "\t" << magy(i,1) << "\t" << magy(i,2) << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==3)//along z
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[2]*geom::dim[2] ; i++)
            {
                sofs << timeid << "\t" << i << "\t" << magz(i,0) << "\t" << magz(i,1) << "\t" << magz(i,2) << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==5)//along x
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[0]*geom::dim[0] ; i++)
            {
                sofs << timeid << "\t" << i;
                for(unsigned int spec = 0 ; spec < geom::ucm.GetNMS() ; spec++)
                {
                    sofs << "\t" << mag_species_x(i,spec,0)*nspl(i,spec) << "\t" << mag_species_x(i,spec,1)*nspl(i,spec) << "\t" << mag_species_x(i,spec,2)*nspl(i,spec);
                }
                sofs << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==6)//along y
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[1]*geom::dim[1] ; i++)
            {
                sofs << timeid << "\t" << i;
                for(unsigned int spec = 0 ; spec < geom::ucm.GetNMS() ; spec++)
                {
                    sofs << "\t" << mag_species_y(i,spec,0)*nspl(i,spec) << "\t" << mag_species_y(i,spec,1)*nspl(i,spec) << "\t" << mag_species_y(i,spec,2)*nspl(i,spec);
                }
                sofs << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==7)//along z
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            //output in index format for plotting with gnuplot
            for(unsigned int i = 0 ; i < geom::Nk[2]*geom::dim[2] ; i++)
            {
                sofs << timeid << "\t" << i;
                for(unsigned int spec = 0 ; spec < geom::ucm.GetNMS() ; spec++)
                {
                    sofs << "\t" << mag_species_z(i,spec,0)*nspl(i,spec) << "\t" << mag_species_z(i,spec,1)*nspl(i,spec) << "\t" << mag_species_z(i,spec,2)*nspl(i,spec);
                }
                sofs << std::endl;
            }
            sofs << std::endl << std::endl;
        }
        else if(spins::mag_calc_method==8)
        {
            const double timeid=static_cast<double>(t)*llg::dt;
            sofs << timeid << "\t" << lx1 << "\t" << ly1 << "\t" << lz1 << "\t" << lx2 << "\t" << ly2 << "\t" << lz2 << std::endl;
        }
        else if(spins::mag_calc_method==9)
        {
            if(disOutForm=="single")
            {
                const double timeid=static_cast<double>(t)*llg::dt;
                for(unsigned int cx = 0 ; cx < maxcx+1 ; cx++)
                {
                    for(unsigned int cy = 0 ; cy < maxcy+1 ; cy++)
                    {
                        for(unsigned int cz = 0 ; cz < maxcz+1 ; cz++)
                        {
                            sofs << timeid << "\t" << cx << "\t" << cy << "\t" << cz << "\t" << magdisc(cx,cy,cz,0) << "\t" << magdisc(cx,cy,cz,1) << "\t" << magdisc(cx,cy,cz,2) << std::endl;
                        }
                    }
                }
            }
            else if(disOutForm=="vtu")
            {
                outputDiscVTU(t);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Method for outputting the discrete magnetization is not recognised. Check flag llg:DiscOutputFormat");
            }
        }

    }


}
