// File: geom.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 16 Jan 2013 16:47:27
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
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
    //the number of atoms in the unit cell
    unsigned int nauc=0;
    //maximum system size
    unsigned int maxss=0;
    //The unit vectors describing the lattice unit cell (unit vectors)
    Array2D<double> L,Linv;
    //The a,b and c values (i.e. lattice constants)
    Array<double> abc;
    //Number of K points
    Array<unsigned int> Nk;
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
        unitCellMembers ucm(nauc);

        config::Info << "Positions:" << std::endl;
        for(unsigned int i = 0 ; i < nauc ; i++)
        {
            std::stringstream sstr;
            sstr << "atom" << i;
            std::string str=sstr.str();
            ucm.SetPosVec(setting[str.c_str()][0],setting[str.c_str()][1],setting[str.c_str()][2],i);
            FIXOUTVEC(config::Info,str,ucm.GetX(i),ucm.GetY(i),ucm.GetZ(i));
        }
        maxss=dim[0]*dim[1]*dim[2]*nauc;
        FIXOUT(config::Info,"Maximum number of spins:" << maxss << std::endl);
        outstruc << maxss << std::endl << std::endl;
        for(unsigned int i = 0 ; i < dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < dim[2] ; k++)
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
                        // the *5 is just a scaling factor
                        outstruc << "Ag\t" << ri[0]*5 << "\t" << ri[1]*5 << "\t" << ri[2]*5 << std::endl;
                    }
                }
            }
        }
        Nk.resize(3);
        abc.resize(3);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            Nk[i]=setting["Nk"][i];
            abc[i]=setting["abc"][i];
        }
        FIXOUT(config::Info,"Number of K-points:" << "[ " << Nk[0] << " , " << Nk[1] << " , " << Nk[2] << " ]" << std::endl);
        FIXOUTVEC(config::Info,"Lattice constants:",abc[0],abc[1],abc[2]);


    }
}
