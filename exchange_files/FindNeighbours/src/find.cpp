// File: find.cpp
// Author:Tom Ostler
// Created: 08 Nov 2015
// Last-modified: 08 Nov 2015 11:37:08

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
    for(unsigned int i = 0 ; i < 3 ; i++)
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
    for(unsigned int i = 0 ; i < 3 ; i++)
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
    for(unsigned int i = 0 ; i < 3 ; i++)
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
    //The last elements store the atom number, the species number
    //and the unit cell that it is part of and whether it is the
    //central cell or not
    coords.resize(nuc[0]*Nk[0],nuc[1]*Nk[1],nuc[2]*Nk[2],4);
    coords.IFill(0);
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
    //lookup returns the atomid, the species and the unit cell that it is in and whether it is in the central cell or not
    lu.resize(ngp,4);
    Array3D<int> uccoords;
    //The unit cell species
    Array3D<int> ucspec;
    Array<int> ucspeclu;
    uccoords.resize(Nk[0],Nk[1],Nk[2]),
    ucspec.resize(Nk[0],Nk[1],Nk[2]);
    ucspeclu.resize(nauc);
    uccoords.IFill(-1);
    ucspec.IFill(-1);
    Info << "================================================" << std::endl;
    for(unsigned int i = 0 ; i < nauc ; i++)
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
        FIXOUTVEC(Info,"Integer coords:",ic[0],ic[1],ic[2]);

        //dump the rest of the information for now
        double dump;
        std::string sdump;
        ucfi >> dump >> dump >> dump;
        ucfi >> sdump;
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
    return(0);
}
