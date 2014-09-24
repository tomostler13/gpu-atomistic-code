// File: geom_glob.cpp
// Author:Tom Ostler
// Created: 26 July 2014
// Last-modified: 24 Sep 2014 12:59:41
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;

//To keep the code tidy the declaration of the geom namespace variables are declared here
namespace geom
{
    //the number of unit cells
    unsigned int dim[3]={0,0,0};
    //zero pad size (in unit cells)
    unsigned int zpdim[3]={0,0,0};
    //For r2c transform the final dimension must be zpdim[2]/2+1
    unsigned int cplxdim=0;
    //number of spins, zeropad size
    //number of magnetic species
    unsigned int nspins=0,zps=0,czps=0,nms=0;
    //The a,b and c values (i.e. lattice constants)
    Array<double> abc,gamma,lambda,llgpf,sublattice,rx,ry,rz;
    //Number of K points
    Array<unsigned int> Nk;
    //lookup array. Give atom number and return coordinates
    Array2D<int> lu;
    //Coords array. Give coords and returns atom number. This coords array
    //is the size of the zero padded arrays and has places where atoms do
    //not exist
    std::string place;
    Array4D<int> coords;
    //instance of the unit cell
    unitCellMembers ucm;

    void readconfig(int argc,char *argv[])
    {
        assert(config::lcf);
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Geometry and magnetic property details***" << std::endl;
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
        FIXOUT(config::Info,"Dipole fields included?" << config::isTF(config::inc_dip) << std::endl);

        libconfig::Setting &setting = config::cfg.lookup("system");
        //do we want to write the unit cell info to the log file?
        bool logunit=false;
        setting.lookupValue("log_unit_cell",logunit);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            //number of unit cells in each direction
            dim[i]=setting["dim"][i];
            zpdim[i]=2*dim[i];
        }
        FIXOUT(config::Info,"Dimensions (unit cells):" << "[ " << dim[0] << " , " << dim[1] << " , " << dim[2] << " ]" << std::endl);
        //string for holding the information about the unit cell file (we do not need global scope here)
        std::string ucf;
        //get the string for opening the unit cell file
        setting.lookupValue("unitcellfile",ucf);
        FIXOUT(config::Info,"Unit cell information file:" << ucf << std::endl);
        //stream for reading unit cell file
        std::ifstream ucfi(ucf.c_str());
        if(!ucfi)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open unit cell file");
        }

        ucfi >> nms;

        FIXOUT(config::Info,"Number of magnetic species (sublattices):" << nms << std::endl);
        unsigned int nauc=0;
        ucfi >> nauc;
        FIXOUT(config::Info,"Number of atoms in the unit cell:" << nauc << std::endl);
        //initialize an instance of the unit cell class
        ucm.init(nauc,nms);
        //loop over the atoms in the unit cell and set the properties in the class (ucm)
        for(unsigned int i = 0 ; i < nauc ; i++)
        {
            unsigned int sl=0;
            //which sublattice do I belong to?
            ucfi >> sl;
            //set the sublattice
            ucm.SetSublattice(sl,i);

            //get the position of this atom in the unit cell
            double c[3]={0,0,0};
            ucfi >> c[0] >> c[1] >> c[2];
            //set the position vector
            ucm.SetPosVec(c[0],c[1],c[2],i);

            double temp=0.0;//this is going to hold mu, lambda and gamma so we call it temp
            //set the magnetic moment
            ucfi >> temp;
            ucm.SetMu(temp,i);
            //set the damping
            ucfi >> temp;
            ucm.SetDamping(temp,i);
            //set the gyromagnetic ratio
            ucfi >> temp;
            ucm.SetGamma(temp,i);
            //Get the element of the atom in the unit cell
            std::string str;
            ucfi >> str;
            ucm.SetElement(str,i);
            //Get the initial spin configuration (reuse the c array)
            ucfi >> c[0] >> c[1] >> c[2];
            ucm.SetInitSpin(c[0],c[1],c[2],i);


        }
        ucfi.close();
        if(ucfi.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close unit cell file.");
        }
        Nk.resize(3);
        abc.resize(3);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            Nk[i]=setting["Nk"][i];
            abc[i]=setting["abc"][i];
        }

    }
}
