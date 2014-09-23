// File: geom.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 23 Sep 2014 13:19:37
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


//All global variables for the geom namespace are in the file
//geom_glob.cpp
namespace geom
{
    void initGeom(int argc,char *argv[])
    {
        assert(config::lcf);
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Geometry details***" << std::endl;
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

        libconfig::Setting &setting = config::cfg.lookup("system");
        std::ofstream outstruc;
        outstruc.open("outstruc.xyz");
        if(!outstruc.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing structure");
        }
        L.resize(3,3);
        L.IFill(0);
        Linv.resize(3,3);
        Linv.IFill(0);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            L(i,i)=1.0;
        }

        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            Linv(0,i)=L(0,i);
            Linv(1,i)=L(1,i);
            Linv(2,i)=L(2,i);
            dim[i]=setting["dim"][i];
            zpdim[i]=2*dim[i];
        }
        //how mamy magnetic types
        FIXOUT(config::Info,"Unit cell matrix (L):" << std::showpos << "[ " << L(0,0) << " , " << L(0,1) << " , " << L(0,2) << " ]" << std::endl);
        FIXOUT(config::Info,"" << "[ " << std::showpos << L(1,0) << " , " << L(1,1) << " , " << L(1,2) << " ]" << std::endl);
        FIXOUT(config::Info,"" << "[ " << std::showpos << L(2,0) << " , " << L(2,1) << " , " << L(2,2) << " ]" << std::endl);
        FIXOUT(config::Info,"Inverting matrix (L):" << std::flush);
        util::inverse(Linv.ptr(),3);
        config::Info << "Done" << std::endl;
        FIXOUT(config::Info,"Unit cell matrix inverse (Linv):" << "[ " << std::showpos << Linv(0,0) << " , " << Linv(0,1) << " , " << Linv(0,2) << " ]" << std::endl);
        FIXOUT(config::Info,"" << "[ " << std::showpos << Linv(1,0) << " , " << Linv(1,1) << " , " << Linv(1,2) << " ]" << std::endl);
        FIXOUT(config::Info,"" << "[ " << std::showpos << Linv(2,0) << " , " << Linv(2,1) << " , " << Linv(2,2) << " ]" << std::endl);
        setting.lookupValue("nauc",nauc);
        FIXOUT(config::Info,"Dimensions (unit cells):" << "[ " << dim[0] << " , " << dim[1] << " , " << dim[2] << " ]" << std::endl);
        FIXOUT(config::Info,"Number of atoms in primitive cell:" << nauc << std::endl);

        ucm.init(nauc);
        config::Info << "Positions:" << std::endl;
        for(unsigned int i = 0 ; i < nauc ; i++)
        {
            std::stringstream sstr;
            sstr << "atom" << i;
            std::string str=sstr.str(),str1;
            ucm.SetPosVec(setting[str.c_str()][0],setting[str.c_str()][1],setting[str.c_str()][2],i);
            FIXOUTVEC(config::Info,str,ucm.GetX(i),ucm.GetY(i),ucm.GetZ(i));

        }
        for(unsigned int i = 0 ; i < nauc ; i++)
        {
            std::stringstream sstr;
            sstr << "Element" << i;
            std::string str=sstr.str(),str1;
            setting.lookupValue(str.c_str(),str1);
            FIXOUT(config::Info,"Element string:" << str1 << std::endl);
            ucm.SetElement(str1,i);
        }

        maxss=dim[0]*dim[1]*dim[2]*nauc;
        Nk.resize(3);
        abc.resize(3);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            Nk[i]=setting["Nk"][i];
            abc[i]=setting["abc"][i];
        }
        FIXOUT(config::Info,"Number of K-points:" << "[ " << Nk[0] << " , " << Nk[1] << " , " << Nk[2] << " ]" << std::endl);
        FIXOUTVEC(config::Info,"Lattice constants:",abc[0],abc[1],abc[2]);
        cplxdim=(dim[2]*Nk[2])+1;

		czps=zpdim[0]*Nk[0]*zpdim[1]*Nk[1]*cplxdim;
		FIXOUT(config::Info,"czps:" << czps << std::endl);
        FIXOUT(config::Info,"Z-dimension for r2c and c2r transforms:" << cplxdim << std::endl);

        nspins=maxss;
        zps=1;
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            zps*=(2*dim[i]*Nk[i]);
        }

        FIXOUT(config::Info,"Maximum number of spins:" << maxss << std::endl);
        FIXOUT(config::Info,"Number of spins:" << nspins << std::endl);
        FIXOUT(config::Info,"Zero pad size:" << zps << std::endl);
        outstruc << maxss << std::endl << std::endl;
        //the 5 entries for each spin correspond to
        // 0,1,2 - the x,y,z positions in the unit cell
        // 3 is the atom number (nothing to do with whether it is magnetic or not) in the unit cell (probably won't be used much
        // 4 is the magnetic species number
        lu.resize(nspins,5);
        lu.IFill(0);
        //the 3 bits of information:
        // 0 - the magnetic atom number
        // 1 - the atom type (not magnetic)
        // 2 - the magnetic species type
        coords.resize(zpdim[0]*Nk[0],zpdim[1]*Nk[1],zpdim[2]*Nk[2],3);
        //-2 here corresponds to empty k-mesh point
        //-1 corresponds to an empty k-mesh point but with an imaginary atom
        //there for the determination of the interaction matrix
        coords.IFill(-2);
        std::ofstream sloc("magat.dat");
        if(sloc.is_open()!=true)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing zero pad information");
        }
        unsigned int atom_counter=0;
        for(unsigned int i = 0 ; i < zpdim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < zpdim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < zpdim[2] ; k++)
                {
                    for(unsigned int t = 0 ; t < nauc ; t++)
                    {
                        double ri[3]={0,0,0};
                        double qi[3]={double(i),double(j),double(k)};
                        for(unsigned int a = 0 ; a < 3 ; a++)
                        {
                            for(unsigned int b = 0 ; b < 3 ; b++)
                            {
                                ri[a]+=L(a,b)*qi[b];
                            }
                            ri[a]+=ucm.GetComp(t,a);
                        }
                        if(i<dim[0] && j<dim[1] && k<dim[2])
                        {
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),0)=atom_counter;
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),1)=t;
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),1)=-1;


                            lu(atom_counter,0)=int(double(Nk[0])*(double(i)+ucm.GetX(t)));
                            lu(atom_counter,1)=int(double(Nk[1])*(double(j)+ucm.GetY(t)));
                            lu(atom_counter,2)=int(double(Nk[2])*(double(k)+ucm.GetZ(t)));
                            lu(atom_counter,3)=t;
                            //set all atoms to zero (i.e initially they are species 0)
                            lu(atom_counter,4)=0;
                            atom_counter++;
                        }
                        else
                        {
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),0)=-1;
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),1)=t;
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),2)=-1;
                        }


                    }
                }
            }
        }
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Magnetic atom placement details***" << std::endl;
        setting.lookupValue("NumberMagneticTypes",nms);
        if(nms > 5)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Too many magnetic types. If you want more edit cuint.cu variable MAXNSPEC (a #define).\n This is because the constant memory allocation for the species dependent variables cannot be\ndynamic.");
        }
        setting.lookupValue("PlaceMagneticType",place);
        FIXOUT(config::Info,"Number of magnetic types:" << nms << std::endl);
        FIXOUT(config::Info,"Method for placing magnetic atoms:" << place << std::endl);
        //check that the atom placement method is recognised
        //random is as it says on the tin
        //NULL is for the case where there is one atom only
        if(place != "random" && place != "NULL" && place!= "unitcell")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Method of placing magnetic atoms not recognised");
        }
        //place each species on the lattice randomly
        if(place=="random")
        {
            setting.lookupValue("
        }

        sloc << "#This file contains the positions of the magnetic species and their type" << std::endl;
        sloc << "# x [m] - y [m] - z[m] - Atom type" << std::endl;
        for(unsigned int i = 0 ; i < geom::dim[0]*Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1]*Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2]*Nk[2] ; k++)
                {
                    if(coords(i,j,k,0)>0)//then the atom exists in this location
                    {
                        sloc << double(i)*abc[0] << "\t" << double(j)*abc[1] << "\t" << double(k)*abc[2] <<
                }
            }
        }
        sloc.close();
        outstruc.close();
        if(sloc.is_open() || outstruc.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close structure files.");
        }
        if(atom_counter!=nspins)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Number of atoms placed on mesh is more or less than expected.");
        }
    }
}
