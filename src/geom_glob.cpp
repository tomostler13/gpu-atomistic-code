// File: geom_glob.cpp
// Author:Tom Ostler
// Created: 26 July 2014
// Last-modified: 08 Sep 2016 16:29:47
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
    //primitive vector
    Array2D<double> rprim;
    //angdeg
    Array<double> angdeg;
    //For r2c transform the final dimension must be zpdim[2]/2+1
    unsigned int cplxdim=0;
    //number of spins, zeropad size
    //number of magnetic species
    unsigned int nspins=0,zps=0,czps=0,nms=0;
    //The a,b and c values (i.e. lattice constants)
    Array<double> abc,mu,gamma,lambda,llgpf,sigma,rx,ry,rz;
    Array<unsigned int> sublattice;
    //Number of mesh points
    Array<unsigned int> Nk;
    //lookup array. Give atom number and return coordinates
    Array2D<int> lu,zplu;
    //Coords array. Give coords and returns atom number. This coords array
    //is the size of the zero padded arrays and has places where atoms do
    //not exist
    std::string place;
    Array4D<int> coords;
    //instance of the unit cell
    unitCellMembers ucm;
    bool logunit=false;
    Array4D<unsigned int> atnolu;

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

        rprim.resize(3,3);rprim.IFill(0);
        //angdeg.resize(3);angdeg.IFill(0);
        libconfig::Setting &setting = config::cfg.lookup("system");
        //do we want to write the unit cell info to the log file?
        setting.lookupValue("log_unit_cell",logunit);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            //number of unit cells in each direction
            dim[i]=setting["dim"][i];
            zpdim[i]=2*dim[i];
        }
        FIXOUT(config::Info,"Dimensions (unit cells):" << "[ " << dim[0] << " , " << dim[1] << " , " << dim[2] << " ]" << std::endl);
        Nk.resize(3);
        abc.resize(3);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            Nk[i]=setting["Nm"][i];
            abc[i]=setting["abc"][i];
        }
        //string for holding the information about the unit cell file (we do not need global scope here)
        std::string ucf;
        //get the string for opening the unit cell file
        setting.lookupValue("unitcellfile",ucf);
        //read rprim
        for(unsigned int i = 0 ; i < 3 ; i++)
        {

            //double magrprim=0.0;
            std::stringstream sstr;
            sstr << "rprim" << i+1;
            std::string str=sstr.str();
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                try
                {
                    rprim(i,j)=setting[str.c_str()][j];
                    rprim(i,j)*=(abc[i]/1e-10);
                    //magrprim+=rprim(i,j)*rprim(i,j);
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Error reading rprim, check and restart");
                }
            }
            FIXOUTVEC(config::Info,str.c_str(),rprim(i,0),rprim(i,1),rprim(i,2));
        }
        /*angdeg[0]=57.8457;
        angdeg[1]=57.8457;
        angdeg[2]=57.8457;
        rprim(0,0)=1.0;
        rprim(0,1)=cos(angdeg[2]*M_PI/180.0);
        rprim(1,1)=sin(angdeg[2]*M_PI/180.0);
        rprim(0,2)=cos(angdeg[1]*M_PI/180.0);
        rprim(1,2)=(cos(angdeg[0]*M_PI/180.0)-rprim(0,1)*rprim(0,2))/rprim(1,1);
        rprim(2,2)=sqrt(1.0-rprim(0,2)*rprim(0,2)-rprim(1,2)*rprim(1,2));

        //rescale
        rprim(0,0)*=(abc[0]/1e-10);
        rprim(0,1)*=(abc[0]/1e-10);
        rprim(0,2)*=(abc[0]/1e-10);
        rprim(1,0)*=(abc[1]/1e-10);
        rprim(1,1)*=(abc[1]/1e-10);
        rprim(1,2)*=(abc[1]/1e-10);
        rprim(2,0)*=(abc[2]/1e-10);
        rprim(2,1)*=(abc[2]/1e-10);
        rprim(2,2)*=(abc[2]/1e-10);

        FIXOUTVEC(config::Info,"rprim1:",rprim(0,0),rprim(0,1),rprim(0,2));
        FIXOUTVEC(config::Info,"rprim2:",rprim(1,0),rprim(1,1),rprim(1,2));
        FIXOUTVEC(config::Info,"rprim3:",rprim(2,0),rprim(2,1),rprim(2,2));*/
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
        //now we know the unit cell size and the dimensions (how many unit cells in each direction) we can allocate atnolu
        atnolu.resize(dim[0],dim[1],dim[2],nauc);
        atnolu.IFill(0);
        unsigned int errcheck=0;//for a bit of error checking from unit cell class
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
            ucm.SetXred(c[0],c[1],c[2],i);
            //set the position vector
            ucm.SetPosVec(static_cast<unsigned int>(c[0]*Nk[0]+0.5),static_cast<unsigned int>(c[1]*Nk[1]+0.5),static_cast<unsigned int>(c[2]*Nk[2]+0.5),i);

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
            //get the first order uniaxial anisotropy constant (reuse the temp variable
            ucfi >> temp;
            //get it's direction and reuse c array
            ucfi >> c[0] >> c[1] >> c[2];
            errcheck=ucm.SetK1U(i,temp,c[0],c[1],c[2]);
            //check what the error is
            if(errcheck!=0)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage(ucm.errorCode(errcheck));
            }


        }
        ucfi.close();
        if(ucfi.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close unit cell file.");
        }

    }
}
