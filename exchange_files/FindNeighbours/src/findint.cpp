// File: find.cpp
// Author:Tom Ostler
// Created: 08 Nov 2015
// Last-modified: 08 Nov 2015 13:48:19

//The purpose of this section of code is to read a unit cell and
//replicate the unit cell out to the specified number in each direction.
//The code will then find the neighbours of each of the atoms in the
//central unit cell and create an exchange file for input into the ASD
//part of the code.
//       |    .     |    .     |    .     |    .     |    .     |
// . . . |    .     |    .     |    .     |    .     |    .     | . . .
//       |    .     |    .     |    .     |    .     |    .     |
//-----------------------------------------------------------------
//       |          |          |          |          |          |
// . . . |          |          |          |          |          | . . .
//       |          |          |          |          |          |
//-----------------------------------------------------------------
//       |          |          |          |          |          |
// . . . |          |          |  Central |          |          | . . .
//       |          |          |    cell  |          |          |
//-----------------------------------------------------------------
//       |          |          |          |          |          |
// . . . |          |          |          |          |          | . . .
//       |          |          |          |          |          |
//-----------------------------------------------------------------
//       |    .     |    .     |    .     |    .     |    .     |
// . . . |    .     |    .     |    .     |    .     |    .     | . . .
//       |    .     |    .     |    .     |    .     |    .     |
#include <cmath>
#include <iostream>
#include <libconfig.h++>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "../../../inc/error.h"
#include "../../../inc/array.h"
#include "../../../inc/array2d.h"
#include "../../../inc/array3d.h"
#include "../../../inc/array4d.h"
#include "../../../inc/defines.h"
int main(int argc,char *argv[])
{
    //for outputting info
    std::ofstream Info("info.out");
    Array<int> Nk;Nk.resize(3);Nk(0)=0;Nk(1)=0;Nk(2)=0;
    Array<int> cd;cd.resize(3);cd(0)=0;cd(1)=0;cd(2)=0;
    Array<double> abc;abc.resize(3);abc(0)=0;abc(1)=0;abc(2)=0;
    if(!Info.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file info.out for outputting information.");
    }
    libconfig::Config cfg;
    if(argc < 2)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("You must give a config file, exiting");
    }
    try
    {
        cfg.readFile(argv[1]);
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
    std::string ucfs;
    libconfig::Setting &setting = cfg.lookup("system");
    if(setting.lookupValue("UnitCellFile",ucfs))
    {
        FIXOUT(Info,"Unit cell file:" << ucfs << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read the unit cell file. Setting system.UnitCellFile (string)");
    }
    //Get the number of grid points
    for(int i = 0 ; i < 3 ; i++)
    {
        try
        {
            Nk[i]=setting["Nk"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
    }
    FIXOUTVEC(Info,"Grid points:",Nk[0],Nk[1],Nk[2]);

    //Get the lattice constants for the distances
    for(int i = 0 ; i < 3 ; i++)
    {
        try
        {
            abc[i]=setting["abc"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
    }
    FIXOUTVEC(Info,"Lattice constants:",abc[0],abc[1],abc[2]);

    //Get the cut-off (in lattice constants)
    for(int i = 0 ; i < 3 ; i++)
    {
        try
        {
            cd[i]=setting["CellDistance"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath();
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }
    }
    FIXOUTVEC(Info,"Lattice cut-off:",cd[0],cd[1],cd[2]);

    //now read the file and put the atoms on the mesh
    //These arrays keep track of the atoms. We are putting
    //the central cell in the middle which requires an off-set

    //Number of unit cells in each direction
    Array<int> nuc;nuc.resize(3);
    nuc[0]=2*cd[0]+1;
    nuc[1]=2*cd[1]+1;
    nuc[2]=2*cd[2]+1;
    Array4D<int> coords;
    Array2D<int> lu;
    Array3D<std::string> atstrings;
    //The last elements store the atom number, the species number
    //and the unit cell that it is part of and whether it is the
    //central cell or not (1 is in cuc and 0 is not)
    coords.resize(nuc[0]*Nk[0],nuc[1]*Nk[1],nuc[2]*Nk[2],4);
    atstrings.resize(nuc[0]*Nk[0],nuc[1]*Nk[1],nuc[2]*Nk[2]);
    coords.IFill(-1);
    std::ifstream ucfi;
    FIXOUT(Info,"Opening unit cell file:" << std::flush);
    ucfi.open(ucfs.c_str());
    if(!ucfi.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open the unit cell file.");
    }
    else
    {
        Info << "Done" << std::endl;
    }
    //number of species
    int nspec=0;
    ucfi >> nspec;
    if(nspec < 1)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Detected zero species, must be at least 1. Check your unit cell file.");
    }
    else
    {
        FIXOUT(Info,"Number of species in the unit cell:" << nspec << std::endl);
    }
    //number of atoms in the unit cell
    int nauc=0;
    ucfi >> nauc;
    if(nauc<1)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Detected no atoms in the unit cell, check your unit cell file.");
    }
    else
    {
        FIXOUT(Info,"Number of atoms in the unit cell:" << nauc << std::endl);
    }

    //the total number of grid points for the lookup
    int ngp=nuc[0]*Nk[0]*nuc[1]*Nk[1]*nuc[2]*Nk[2];
    FIXOUT(Info,"Total number of grid points:" << ngp << std::endl);
    //lookup returns the global atom coords, the species and the unit cell that it is in and whether it is in the central cell or not
    lu.resize(ngp,6);
    Array3D<int> uccoords;
    //The unit cell species
    Array3D<int> ucspec;
    Array<int> ucspeclu;
    Array2D<int> uclu;
    Array<std::string> ucstrlu;
    ucstrlu.resize(nauc);
    uclu.resize(nauc,3);
    uclu.IFill(0);
    uccoords.resize(Nk[0],Nk[1],Nk[2]),
        ucspec.resize(Nk[0],Nk[1],Nk[2]);
    ucspeclu.resize(nauc);
    uccoords.IFill(-1);
    ucspec.IFill(-1);
    Info << "================================================" << std::endl;
    for(int i = 0 ; i < nauc ; i++)
    {
        ucfi >> ucspeclu[i];
        double lc[3]={0,0,0};
        ucfi >> lc[0] >> lc[1] >> lc[2];
        std::stringstream sstr;
        sstr << "Atom number " << i << " is at location:";
        std::string str=sstr.str();
        FIXOUTVEC(Info,str,lc[0],lc[1],lc[2]);
        int ic[3]={0,0,0};
        ic[0]=static_cast<int>(lc[0]*static_cast<double>(Nk[0])+0.5);
        ic[1]=static_cast<int>(lc[1]*static_cast<double>(Nk[1])+0.5);
        ic[2]=static_cast<int>(lc[2]*static_cast<double>(Nk[2])+0.5);
        uccoords(ic[0],ic[1],ic[2])=i;
        uclu(i,0)=ic[0];
        uclu(i,1)=ic[1];
        uclu(i,2)=ic[2];
        FIXOUTVEC(Info,"Integer coords:",ic[0],ic[1],ic[2]);

        //dump the rest of the information for now
        double dump;
        std::string sdump;
        ucfi >> dump >> dump >> dump;
        ucfi >> ucstrlu[i];
        ucfi >> dump >> dump >> dump >> dump >> dump >> dump >> dump;
    }
    Info << "================================================" << std::endl;
    //work out which is the central unit cell
    int cuc[3]={0,0,0};
    cuc[0]=cd[0];
    cuc[1]=cd[1];
    cuc[2]=cd[2];
    FIXOUTVEC(Info,"Central unit cell coords:",cuc[0]+1,cuc[1]+1,cuc[2]+1);
    FIXOUTVEC(Info,"(C array) central unit cell:",cuc[0],cuc[1],cuc[2]);


    //loop over the unit cells and place the atoms on the grid
    int atom_counter=0;
    //unit cell count
    int uccount=0;
    //debug - check the number of atoms in the central cell. Should be equal to nauc
    int check=0;
    //debug - check that an atom is not placed at the same mesh coordinate (duplicated atom)
    Array3D<int> checkcoord;
    checkcoord.resize(nuc[0]*Nk[0],nuc[1]*Nk[1],nuc[2]*Nk[2]);
    checkcoord.IFill(0);
    for(int i = 0 ; i < nuc[0] ; i++)
    {
        for(int j = 0 ; j < nuc[1] ; j++)
        {
            for(int k = 0 ; k < nuc[2] ; k++)
            {
                for(int aiuc = 0 ; aiuc < nauc ; aiuc++)
                {
                    int lc[3]={0,0,0};
                    lc[0]=i*Nk[0]+uclu(aiuc,0);
                    lc[1]=j*Nk[1]+uclu(aiuc,1);
                    lc[2]=k*Nk[2]+uclu(aiuc,2);
                    if(checkcoord(lc[0],lc[1],lc[2])==1)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream sstr;
                        //std::cerr << "\n" << atom_counter << "\t" << i << "\t" << j << "\t" << k << "\t" << lc[0] << "\t" << lc[1] << "\t" << lc[2] <<std::endl;
                        //std::cerr << uclu(aiuc,0) << "\t" << uclu(aiuc,1) << "\t" << uclu(aiuc,2) << std::endl;
                        //std::cerr << aiuc << "\t" << coords(lc[0],lc[1],lc[2],0) << std::endl;
                        sstr << "A duplicate atom has been detected with coords (integer):\t" << i << " , " << j << " , " << k << std::endl;
                        std::string str=sstr.str();
                        error::errMessage(str);
                    }

                    coords(lc[0],lc[1],lc[2],0)=atom_counter;
                    coords(lc[0],lc[1],lc[2],1)=ucspeclu[aiuc];
                    coords(lc[0],lc[1],lc[2],2)=uccount;
                    lu(atom_counter,0)=lc[0];
                    lu(atom_counter,1)=lc[1];
                    lu(atom_counter,2)=lc[2];
                    lu(atom_counter,3)=ucspeclu[aiuc];
                    lu(atom_counter,4)=uccount;
                    atstrings(lc[0],lc[1],lc[2])=ucstrlu(aiuc);
                    checkcoord(i,j,k)=true;
                    //check if we are in the central unit cell
                    if(i==cuc[0] && j==cuc[1] && k==cuc[2])
                    {
                        coords(lc[0],lc[1],lc[2],3)=1;
                        lu(atom_counter,5)=1;
                        check++;
                    }
                    checkcoord(lc[0],lc[1],lc[2])=1;
                    atom_counter++;
                }

                uccount++;
            }
        }
    }
    if(check!=nauc)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The number of atoms in the central unit cell is not consistent");
    }

    //loop over all of the atoms in the (central) unit cell
    for(int aiuc = 0 ; aiuc < nauc ; aiuc++)
    {
        std::stringstream sstr;
        sstr << "atom" << aiuc << "_interactions.dat";
        std::string str=sstr.str();
        std::ofstream ofs(str.c_str());
        if(!ofs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Could not open file for atom " << aiuc << " in the unit cell.";
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }

        ofs << "#Atom i id - atom type (label) - neighbour j id - neigh atom type (label) - species i - species j - unit cell of i - unit cell of j - {integer coords of i} - {integer coords of j} - {real space interaction vector} - distance - {mesh interaction vector}" << std::endl;
        //the current atom in the unit cell's position
        int lc[3]={0,0,0};
        lc[0]=cuc[0]*Nk[0]+uclu(aiuc,0);
        lc[1]=cuc[1]*Nk[1]+uclu(aiuc,1);
        lc[2]=cuc[2]*Nk[2]+uclu(aiuc,2);
        //Get the atom number
        int atid=coords(lc[0],lc[1],lc[2],0);
        //check that indeed this atom is in the central unit cell.
        //This is a little redundant but good for checking anyway
        if(coords(lc[0],lc[1],lc[2],3)==0)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("This is a bug. The code is trying to look for the neighbours of an atom that is not in the central unit cell (this should not be allowed to happen).");
        }
        //loop over all of the atoms in the system up to the distance specified
        //in the input file. Loop over space rather than atoms (less efficient
        //this way but oh well).
        for(int i = 0 ; i < nuc[0]*Nk[0] ; i++)
        {
            for(int j = 0 ; j < nuc[1]*Nk[1] ; j++)
            {
                for(int k = 0 ; k < nuc[2]*Nk[2] ; k++)
                {
                    //check first if the atom exists
                    if(coords(i,j,k,0)>-1)
                    {
                        //another check that it is not the same atom (i.e. no self interaction)
                        if(i==lc[0] && j==lc[1] && k==lc[2])
                        {
                            //then it is a self interaction and we should skip it
                        }
                        else
                        {
                            //calculate the lookup (in k-points)
                            int luvec[3]={lc[0]-i,lc[1]-j,lc[2]-k};
                            //convert that to real space
                            double dluvec[3]={0,0,0};
                            for(int xyz = 0 ; xyz < 3 ; xyz++)
                            {
                                dluvec[xyz]=abc[xyz]*static_cast<double>(luvec[xyz])/(static_cast<double>(Nk[xyz]));
                            }//end of xyz loop
                            //neighbour id
                            int nid = coords(i,j,k,0);
                            ofs << atid << "\t" << atstrings(lc[0],lc[1],lc[2]) << "\t" << nid << "\t" << atstrings(i,j,k) << "\t" << coords(lc[0],lc[1],lc[2],1) << "\t" << coords(i,j,k,1) << "\t" << coords(lc[0],lc[1],lc[2],2) << "\t" << coords(i,j,k,2) << "\t[ " << dluvec[0] << " , " << dluvec[1] << " , " << dluvec[2] << " ] " << sqrt(dluvec[0]*dluvec[0]+dluvec[1]*dluvec[1]+dluvec[2]*dluvec[2]) << "\t[ " << luvec[0] << " , " << luvec[1] << " , " << luvec[2] << " ]" << std::endl;
                        }//end of self interaction check
                    }//end of checking if atom exists check
                }
            }
        }
        ofs.close();
        if(ofs.is_open())
        {
            std::stringstream errsstr;
            errsstr << "Could not close output file for atom " << aiuc << " in the unit cell.";
            std::string errstr=errsstr.str();
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage(errstr);
        }
//        exit(0);
    }//aiuc
    return(0);
}
