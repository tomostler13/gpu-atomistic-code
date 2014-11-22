// File: make_chess_board.cpp
// Author:Tom Ostler
// Created: 22 Nov 2014
// Last-modified: 22 Nov 2014 18:12:25

//The purpose of this section of code is to create a unit cell
//file for use with the main program. The specific type of unit
//cell that this code creates is a disordered system of two
//species (the initial idea was to create Gd/Fe). Furthermore,
//the chess board idea comes from creating areas with different
//concentrations. At the moment there are only two concentrations.
//The code could be generalized to N but at the moment it isn't
//required
//       |    .     |    .     |    .     |    .     |    .     |
// . . . |    .     |    .     |    .     |    .     |    .     | . . .
//       |    .     |    .     |    .     |    .     |    .     |
//-----------------------------------------------------------------
//       |          |          |          |          |          |
// . . . |  20% Gd  |  30% Gd  |  20% Gd  |  30% Gd  |  20 %Gd  | . . .
//       |          |          |          |          |          |
//-----------------------------------------------------------------
//       |          |          |          |          |          |
// . . . |  30% Gd  |  20% Gd  |  30% Gd  |  20% Gd  |  30 %Gd  | . . .
//       |          |          |          |          |          |
//-----------------------------------------------------------------
//       |          |          |          |          |          |
// . . . |  20% Gd  |  30% Gd  |  20% Gd  |  30% Gd  |  20 %Gd  | . . .
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
    std::ofstream Info("structure_info.dat");
    unsigned int dim[3]={0,0,0},globdim[3]={0,0,0},nk[3]={1,1,1},nms = 1,nauc=1;
    Array<double> damp,gamma,anis;
    Array2D<double> basis,anisvec,spin;
    Array2D<double> amounts;
    std::vector<std::string> spec;

    if(!Info.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file structure_info.dat for outputting information.");
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
    libconfig::Setting &setting = cfg.lookup("system");
    bool errstatus=false;
    errstatus=setting.lookupValue("nauc",nauc);
    if(errstatus)
    {
        FIXOUT(Info,"Number of atoms in the unit cell (ignoring magnetic atoms):" << nauc << std::endl);
        basis.resize(nauc,3);
        basis.IFill(0);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read number of atoms in the unit cell.");
    }
    errstatus=setting.lookupValue("nms",nms);

    Info << "---------------------------------------------------------------" << std::endl;
    Info << "            Magnetic Species Info                              " << std::endl;
    if(errstatus)
    {
        FIXOUT(Info,"Number of magnetic species:" << nms << std::endl);

        amounts.resize(2,nms);amounts.IFill(0);
        damp.resize(nms);damp.IFill(0);
        gamma.resize(nms);gamma.IFill(0);
        anis.resize(nms);anis.IFill(0);
        anisvec.resize(nms,3);anisvec.IFill(0);
        spec.resize(nms);//leave uninitialized!!!
        spin.resize(nms,3);spin.IFill(0);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read number of atoms in the unit cell.");
    }
    //loop over the magnetic species
    for(unsigned int i = 0 ; i < nms ; i++)
    {
        Info << "---------------------------------------------------------------" << std::endl;
        FIXOUT(Info,"Species:" << i << std::endl);
        damp[i]=setting["damping"][i];
        gamma[i]=setting["gamma"][i];
        anis[i]=setting["D0"][i];
        std::string temp=setting["spec"][i];
        spec[i]=temp;
        //the first index here is the first area
        amounts(0,i)=setting["amount0"][i];
        amounts(1,i)=setting["amount1"][i];
        std::stringstream sstr1;sstr1 << "Amount of species " << i << " in area 0:";
        std::string str1=sstr1.str();
        FIXOUT(Info,str1 << amounts(0,i) << " %" << std::endl);
        sstr1.str("");sstr1 << "Amount of species " << i << " in area 1:";str1=sstr1.str();
        FIXOUT(Info,str1 << amounts(1,i) << " %" << std::endl);
        FIXOUT(Info,"Damping:" << damp[i] << std::endl);
        FIXOUT(Info,"Gamma:" << gamma[i] << std::endl);
        FIXOUT(Info,"Anisotropy:" << anis[i] << std::endl);
        FIXOUT(Info,"String for magnetic atom:" << spec[i] << std::endl);

        std::stringstream sstr;sstr << "spin" << i;std::string str=sstr.str();
        spin(i,0)=setting[str.c_str()][0];
        spin(i,1)=setting[str.c_str()][1];
        spin(i,2)=setting[str.c_str()][2];
        FIXOUTVEC(Info,"Initial spin position:",spin(i,0),spin(i,1),spin(i,2));
        sstr.str("");sstr << "DVec" << i;str=sstr.str();
        anisvec(i,0)=setting[str.c_str()][0];
        anisvec(i,1)=setting[str.c_str()][1];
        anisvec(i,2)=setting[str.c_str()][2];
        FIXOUTVEC(Info,"Anisotropy axis:",anisvec(i,0),anisvec(i,1),anisvec(i,2));
    }
    if(fabs(100-(amounts(0,0)+amounts(0,1)))>1e-12)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The amounts in area 0 do not equal 100%");
    }
    if(fabs(100-(amounts(1,0)+amounts(1,1)))>1e-12)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The amounts in area 1 do not equal 100%");
    }
    Info << "---------------------------------------------------------------" << std::endl;
    Info << "            Atomic/Structural Info                             " << std::endl;
    //loop over the unit cells and read the basis vectors
    for(unsigned int i = 0 ; i < nauc ; i++)
    {
        std::stringstream sstr;
        sstr << "atom" << i;
        std::string str=sstr.str();
        //loop over xyz coords
        for(unsigned int j = 0 ; j < 3 ; j++)
        {
            basis(i,j)=setting[str.c_str()][j];
        }
        FIXOUTVEC(Info,"Basis vectors:",basis(i,0),basis(i,1),basis(i,2));
    }
    for(unsigned int i = 0 ; i < 3 ; i++)
    {
        dim[i]=setting["dim"][i];
        globdim[i]=setting["globdim"][i];
        nk[i]=setting["nk"][i];
    }
    FIXOUTVEC(Info,"Number of unit cells per chess board volume:",dim[0],dim[1],dim[2]);
    FIXOUTVEC(Info,"Number of chess board volumes:",globdim[0],globdim[1],globdim[2]);
    FIXOUTVEC(Info,"Number of mesh points per unit cell:",nk[0],nk[1],nk[2]);

    Array3D<int> coords;
    //arguements here are kx,ky,kz and the return value is the species number (-1 means no atom)
    coords.resize(dim[0]*globdim[0]*nk[0],dim[1]*globdim[1]*nk[1],dim[2]*globdim[2]*nk[2]);
    coords.IFill(-1);
    unsigned long int atom_counter=0;
    //loop over every unit cell and make all atoms (that exist) species 0
    for(unsigned int i = 0 ; i < globdim[0]*dim[0] ; i++)
    {
        for(unsigned int j = 0 ; j < globdim[1]*dim[1] ; j++)
        {
            for(unsigned int k = 0 ; k < globdim[2]*dim[2] ; k++)
            {
                for(unsigned int ca = 0 ; ca < nauc ; ca++)
                {
                    coords(i*nk[0]+static_cast<unsigned int>(static_cast<double>(nk[0])*basis(ca,0)+0.5),
                           j*nk[1]+static_cast<unsigned int>(static_cast<double>(nk[1])*basis(ca,1)+0.5),
                           k*nk[2]+static_cast<unsigned int>(static_cast<double>(nk[2])*basis(ca,2)+0.5))=0;
                    atom_counter++;

                }
            }
        }
    }
    std::ofstream atstruc("atom_struc.xyz");
    if(!atstruc.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open (raw) atomic structure file.");
    }
    atstruc << atom_counter << std::endl << std::endl;
    for(unsigned int i = 0 ; i < globdim[0]*dim[0]*nk[0] ; i++)
    {
        for(unsigned int j = 0 ; j < globdim[1]*dim[1]*nk[1] ; j++)
        {
            for(unsigned int k = 0 ; k < globdim[2]*dim[2]*nk[2] ; k++)
            {
                if(coords(i,j,k)==0)
                {
                    atstruc << "H" << "\t" << 3*i << "\t" << 3*j << "\t" << 3*k << std::endl;
                }
            }
        }
    }
    atstruc.close();
    if(atstruc.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not close the (raw) atom structure file.");
    }
    return(EXIT_SUCCESS);
}
