// File: gen_dw_inp.cpp
// Author:Tom Ostler
// Last-modified: 19 Aug 2022 11:58:26
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <libconfig.h++>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <string>
#include <cstdarg>
#include "../inc/error.h"
#include "../../inc/defines.h"
#include "../../inc/array.h"
#include "../../inc/array2d.h"
#include "../../inc/array3d.h"
#include "../../inc/array4d.h"
#include "../../inc/unitcell.h"
int main(int argc,char *argv[])
{
    if(argc < 7)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("You must give a config file followed by width - centre - direction of DW (0=x, 1=y,2=z) - direction of rotation of spins (0=x, 1=y,2=z) - init spin alignment (0=x, 1=y,2=z) , exiting");
    }

    double delta=atof(argv[2]);
    double width=atof(argv[3]);
    int dirwall=atoi(argv[4]);
    int dirrot=atoi(argv[5]);
    int initspinalign=atoi(argv[6]);
    FIXOUT(std::cout,"DW initial width (delta):" << delta << " [A]" << std::endl);
    FIXOUT(std::cout,"DW initial position (b):" << delta << " [A]" << std::endl);
    FIXOUT(std::cout,"Direction of domain wall:" << dirwall << std::endl);
    FIXOUT(std::cout,"Direction of rotation of domain wall:" << dirrot << std::endl);
    FIXOUT(std::cout,"Initial spin alignment:" << initspinalign << std::endl);


    libconfig::Config cfg;
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
        //exit(EXIT_FAILURE);
    }
    catch(...)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Unspecified libconfig error on reading the input file.");
    }

    Array<int> Nk;
    Array<int> dim;
    Array<double> abc;
    Nk.resize(3);
    Nk.IFill(0);
    abc.resize(3);
    abc.IFill(0.0);
    dim.resize(3);
    dim.IFill(0);

    //READ THE GEOMETRY DETAILS FROM THE CONFIG FILE
    libconfig::Setting &sysset = cfg.lookup("system");
    for(unsigned int i = 0 ; i < 3 ; i++)
    {

        //number of unit cells in each direction
        try
        {
            dim[i]=sysset["dim"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " must be set.";
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }

        int count=3;
        try
        {
            Nk[i]=sysset["Nm"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " element " << i;
            std::string errstr=errsstr.str();
            error::errWarning(errstr);
        }
        try
        {
            abc[i]=sysset["abc"][i];
        }
        catch(const libconfig::SettingNotFoundException &snf)
        {
            error::errPreamble(__FILE__,__LINE__);
            std::stringstream errsstr;
            errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " must be set.";
            std::string errstr=errsstr.str();
            error::errMessage(errstr);
        }

    }
    FIXOUTVEC(std::cout,"Nk = ",Nk[0],Nk[1],Nk[2]);
    FIXOUTVEC(std::cout,"abc = ",abc[0],abc[1],abc[2]);
    FIXOUTVEC(std::cout,"Dimensions (unit cells):",dim[0], dim[1], dim[2]);
    //string for holding the information about the unit cell file (we do not need global scope here)
    std::string ucf;

    Array2D<double> rprim;
    rprim.resize(3,3);
    rprim.IFill(0.0);
    bool rprimset=false;

    //read rprim
    for(unsigned int i = 0 ; i < 3 ; i++)
    {
        int count=9;
        //double magrprim=0.0;
        std::stringstream sstr;
        sstr << "rprim" << i+1;
        std::string str=sstr.str();
        for(unsigned int j = 0 ; j < 3 ; j++)
        {
            try
            {
                rprim(i,j)=sysset[str.c_str()][j];
                rprim(i,j)*=(abc[i]/1e-10);
                //magrprim+=rprim(i,j)*rprim(i,j);
            }
            catch(const libconfig::SettingNotFoundException &snf)
            {
                count--;
                error::errPreamble(__FILE__,__LINE__);
                std::stringstream errsstr;
                errsstr << "Setting not found exception caught. Setting " << snf.getPath();
                std::string errstr=errsstr.str();
                error::errWarning(errstr);
            }
        }
        if(count<9)
        {
            rprimset=false;
        }
        else
        {
            rprimset=true;
        }
        FIXOUTVEC(std::cout,str.c_str(),rprim(i,0),rprim(i,1),rprim(i,2));
    }


    //get the string for opening the unit cell file
    if(!sysset.lookupValue("Unitcellfile",ucf))
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("A unit cell file must be specified (system.Unitcellfile (string))");
    }

    std::ifstream ucfi(ucf.c_str());
    if(!ucfi)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open unit cell file");
    }
    int nms=0,nauc=0;
    unitCellMembers ucm;
    FIXOUT(std::cout,"Unit cell file:" << ucf << std::endl);
    ucfi >> nms;
    FIXOUT(std::cout,"Number of magnetic species (sublattices):" << nms << std::endl);
    nauc=0;
    ucfi >> nauc;
    if(nauc<1)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("There are no atoms in the unit cell, check your unit cell file.");
    }
    FIXOUT(std::cout,"Number of atoms in the unit cell:" << nauc << std::endl);
    //initialize an instance of the unit cell class
    ucm.init(nauc,nms);
    Array4D<unsigned int> atnolu;
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
        //If Nkset is true then we can put the moments on a mesh
        //which is specified in the unit cell file
        int mp[3]={0,0,0};
        ucfi >> c[0] >> c[1] >> c[2];
        ucfi >> mp[0] >> mp[1] >> mp[2];
        ucm.SetMP(mp[0],mp[1],mp[2],i);
        ucm.SetXred(c[0],c[1],c[2],i);

        double temp=0.0;//this is going to hold mu, lambda and gamma so we call it temp
        //set the magnetic moment
        ucfi >> temp;
        ucm.SetMu(temp,i);
        if(temp<1e-5)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Magnetic moment is less than 1e-5 (so, very small). Check your unit cell file is of the correct format");
        }
        else if(temp>1e5)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Magnetic moment is greater than 1e5 (so, really huge). Check your unit cell file is of the correct format");
        }
        //set the damping
        ucfi >> temp;
        ucm.SetDamping(temp,i);
        //set the gyromagnetic ratio
        ucfi >> temp;
        ucm.SetGamma(temp,i);
        if(temp<1e-5)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Gyromagnetic is less than 1e-5 (so, really small). Check your unit cell file is of the correct format");
        }
        else if(temp>1e5)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Gyromagnetic is greater than 1e5 (so, really huge). Remember it should be in units of the free electron value 1.76e11.\nCheck your unit cell file is of the correct format");
        }
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

    std::ofstream ofs("dw_grid.in");
    if(!ofs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open file for outputting spin config");
    }
    int nspins=dim[0]*dim[1]*dim[2]*ucm.NumAtomsUnitCell();
    FIXOUT(std::cout,"Number of spins:" << nspins << std::endl);
    ofs << nspins << std::endl;

    Array<double> init_sign;
    init_sign.resize(nauc);
    init_sign.IFill(0.0);
    for(unsigned int t = 0 ; t < ucm.NumAtomsUnitCell() ; t++)
    {
        init_sign(t)=ucm.GetInitS(t,initspinalign);
    }

    for(unsigned int i = 0 ; i < dim[0] ; i++)
    {
        for(unsigned int j = 0 ; j < dim[1] ; j++)
        {
            for(unsigned int k = 0 ; k < dim[2] ; k++)
            {
                for(unsigned int t = 0 ; t < ucm.NumAtomsUnitCell() ; t++)
                {
                    int gp[3]={0,0,0};
                    gp[0]=i*Nk[0]+ucm.GetMP(t,0);
                    gp[1]=j*Nk[1]+ucm.GetMP(t,1);
                    gp[2]=k*Nk[2]+ucm.GetMP(t,2);
                    ofs << gp[0] << "\t" << gp[1] << "\t" << gp[2] << "\t";

                    // The real space position is equal to
                    // r_x = ((k_x+i*N_{k,x})/N_{k,x})*lattice constant_x
                    // where k_x is the location of the atom in the unit cell
                    // N_{k,x} is the number of mesh-points in the x direction
                    // This can of course be generalised the each direction
                    double xred[3]={ucm.GetXred(t,0),ucm.GetXred(t,1),ucm.GetXred(t,2)};

                    double xcart[3]={0,0,0};
                    xcart[0]=(xred[0]+static_cast<double>(i))*rprim(0,0)+(xred[1]+static_cast<double>(j))*rprim(1,0)+(xred[2]+static_cast<double>(k))*rprim(2,0);
                    xcart[1]=(xred[0]+static_cast<double>(i))*rprim(0,1)+(xred[1]+static_cast<double>(j))*rprim(1,1)+(xred[2]+static_cast<double>(k))*rprim(2,1);
                    xcart[2]=(xred[0]+static_cast<double>(i))*rprim(0,2)+(xred[1]+static_cast<double>(j))*rprim(1,2)+(xred[2]+static_cast<double>(k))*rprim(2,2);

                    double s[3]={0.0,0.0,0.0};
                    s[initspinalign]=-init_sign(t)*tanh((xcart[dirwall]-width)/delta);
                    s[dirrot]=init_sign(t)*sqrt(1.0-s[initspinalign]*s[initspinalign]);
                    ofs << s[0] << "\t" << s[1] << "\t" << s[2] << std::endl;
                }
            }
        }
    }
    ofs.close();
    if(ofs.is_open())
    {
        error::errWarnPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close output file dw_grid.in");
    }




    return(EXIT_SUCCESS);
}
