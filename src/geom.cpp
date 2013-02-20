// File: geom.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 20 Feb 2013 12:44:34
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


namespace geom
{

    //the number of unit cells
    unsigned int dim[3]={0,0,0};
    //zero pad size (in unit cells)
    unsigned int zpdim[3]={0,0,0};
    //For r2c transform the final dimension must be zpdim[2]/2+1
    unsigned int cplxdim=0;
    //the number of atoms in the unit cell
    unsigned int nauc=0;
    //maximum system size
    unsigned int maxss=0,nspins=0,zps=0,czps=0;
    //The unit vectors describing the lattice unit cell (unit vectors)
    Array2D<double> L,Linv;
    //The a,b and c values (i.e. lattice constants)
    Array<double> abc;
    //Number of K points
    Array<unsigned int> Nk;
    //lookup array. Give atom number and return coordinates
    Array2D<int> lu;
    //Coords array. Give coords and returns atom number. This coords array
    //is the size of the zero padded arrays and has places where atoms do
    //not exist
    Array4D<int> coords;
    //instance of the unit cell
    unitCellMembers ucm;
    //check if we are zeropadding the fourier transform
    bool zpcheck=true;
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
        setting.lookupValue("zpcheck",zpcheck);
        if(zpcheck==false && config::incdip==true)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You must zero pad the fourier transform if you want to include the dipolar field terms.");
        }
        FIXOUT(config::Info,"Zero pad the fourier transform:" << config::isTF(zpcheck) << std::endl);
        L.resize(3,3);
        L.IFill(0);
        Linv.resize(3,3);
        Linv.IFill(0);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            L(0,i)=setting["Lex"][i];
            L(1,i)=setting["Ley"][i];
            L(2,i)=setting["Lez"][i];
            Linv(0,i)=L(0,i);
            Linv(1,i)=L(1,i);
            Linv(2,i)=L(2,i);
        }
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            dim[i]=setting["dim"][i];
            if(zpcheck==true)
            {
                zpdim[i]=2*dim[i];
            }
            else
            {
                zpdim[i]=dim[i];
            }
        }
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
        cplxdim=(zpdim[2]*Nk[2]/2)+1;

		czps=zpdim[0]*Nk[0]*zpdim[1]*Nk[1]*cplxdim;
		FIXOUT(config::Info,"czps:" << czps << std::endl);
        FIXOUT(config::Info,"Z-dimension for r2c and c2r transforms:" << cplxdim << std::endl);

        nspins=maxss;
        zps=1;
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            zps*=(zpdim[i]*Nk[i]);
        }

        FIXOUT(config::Info,"Maximum number of spins:" << maxss << std::endl);
        FIXOUT(config::Info,"Number of spins:" << nspins << std::endl);
        FIXOUT(config::Info,"Zero pad size:" << zps << std::endl);
        outstruc << maxss << std::endl << std::endl;
        lu.resize(nspins,4);
        lu.IFill(0);
        coords.resize(zpdim[0]*Nk[0],zpdim[1]*Nk[1],zpdim[2]*Nk[2],2);
        //-2 here corresponds to empty k-mesh point
        //-1 corresponds to an empty k-mesh point but with an imaginary atom
        //there for the determination of the interaction matrix
        coords.IFill(-2);
        std::ofstream sloc("zeropad.dat");
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

                            lu(atom_counter,0)=int(double(Nk[0])*(double(i)+ucm.GetX(t)));
                            lu(atom_counter,1)=int(double(Nk[1])*(double(j)+ucm.GetY(t)));
                            lu(atom_counter,2)=int(double(Nk[2])*(double(k)+ucm.GetZ(t)));
                            lu(atom_counter,3)=t;
                            atom_counter++;
                        }
                        else
                        {
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),0)=-1;
                            coords(int(double(Nk[0])*(double(i)+ucm.GetX(t))),int(double(Nk[1])*(double(j)+ucm.GetY(t))),int(double(Nk[2])*(double(k)+ucm.GetZ(t))),1)=t;
                        }


                    }
                }
            }
        }
        for(unsigned int i = 0 ; i < geom::zpdim[0]*Nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::zpdim[1]*Nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::zpdim[2]*Nk[2] ; k++)
                {
                    if(coords(i,j,k,0)>=0)
                    {
                        sloc << "\"" << ucm.GetElement(coords(i,j,k,1)) << "\"\t" << i << "\t" << j << "\t" << k << std::endl;
                        // the *5 is just a scaling factor
                        outstruc << ucm.GetElement(coords(i,j,k,1)) << "\t" << i*2 << "\t" << j*2 << "\t" << k*2 << std::endl;
                    }
                    else
                    {
                        sloc << "\".\"\t" << i << "\t" << j << "\t" << k << std::endl;
                    }
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
