// File: geom.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 06 Oct 2014 19:21:43
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/geom.h"
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/unitcell.h"
#include "../inc/anis.h"
#include "../inc/llg.h"
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
        readconfig(argc,argv);
        //The unit cell file has a strict format
        //first line -  number of distinct atomic type
        //second line - number of magnetic atoms in unit cell
        //sublattice - kx - ky - kz - mu - lambda - gamma - Element - sx - sy - sz

        //calculate the number of spins
        nspins=dim[0]*dim[1]*dim[2]*ucm.NumAtomsUnitCell();
        //resize these 1D arrays. The atom number should return the value
        mu.resize(nspins);gamma.resize(nspins);lambda.resize(nspins);sigma.resize(nspins);llgpf.resize(nspins);rx.resize(nspins);ry.resize(nspins);rz.resize(nspins);sublattice.resize(nspins);
        //resize the anisotropy arrays
        anis::k1u.resize(nspins);anis::k1udir.resize(nspins,3);
        FIXOUTVEC(config::Info,"Number of K-points:",Nk[0],Nk[1],Nk[2]);
        FIXOUTVEC(config::Info,"Lattice constants:",abc[0],abc[1],abc[2]);
        //For real to complex (or c2r) transforms we can save computation
        //by exploiting half dim size of this type of transform
        cplxdim=(dim[2]*Nk[2])+1;

        //The total number of elements involved in the 3D transform
        czps=zpdim[0]*Nk[0]*zpdim[1]*Nk[1]*cplxdim;
        FIXOUT(config::Info,"czps:" << czps << std::endl);
        FIXOUT(config::Info,"Z-dimension for r2c and c2r transforms:" << cplxdim << std::endl);

        //Calculate the size of the real space zero padded array
        zps=1;
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            zps*=(2*dim[i]*Nk[i]);
        }
        //The point of this array is to number each k-point and assign it a value on the k-point
        //mesh. The reason is so that when we do the convolution on the GPU we can look up the
        //correct array elements
        zplu.resize(zps,3);
        FIXOUT(config::Info,"Number of spins:" << nspins << std::endl);
        FIXOUT(config::Info,"Zero pad size:" << zps << std::endl);
        //Next we want to output the information stored in the class to an output file.
        //If the number of atoms in the unit cell is large (usually when we are using a bit supercell)
        //the we don't want to mess up the format of the output file as it would make it unreadable if
        //we output the information for potentially millions of spins.
        //We will always output the information for the first 5 atoms but the whole data will
        //be output to the log file if selected in the config file.
        if(ucm.NumAtomsUnitCell() > 5)
        {
            config::openLogFile();
        }
        config::printline(config::Info);
        for(unsigned int i = 0 ; i < ucm.NumAtomsUnitCell() ; i++)
        {
            //output to the config file (up to 5)
            if(i < 5)
            {
                FIXOUT(config::Info,"Unit cell atom:" << i << std::endl);
                FIXOUT(config::Info,"Element:" << ucm.GetElement(i) << std::endl);
                FIXOUTVEC(config::Info,"Positions [real space]:",ucm.GetCoord(i,0)/static_cast<double>(Nk[0]),ucm.GetCoord(i,1)/static_cast<double>(Nk[1]),ucm.GetCoord(i,2)/static_cast<double>(Nk[2]));
                FIXOUTVEC(config::Info,"Positions [k-lattice]:",ucm.GetCoord(i,0),ucm.GetCoord(i,1),ucm.GetCoord(i,2));
                FIXOUT(config::Info,"Part of sublattice:" << ucm.GetSublattice(i) << std::endl);
                FIXOUT(config::Info,"Damping:" << ucm.GetDamping(i) << std::endl);
                FIXOUT(config::Info,"Gyromagnetic ratio:" << ucm.GetGamma(i) << " [gamma_free]" << std::endl);
                FIXOUT(config::Info,"Magnetic moment:" << ucm.GetMu(i) << " [muB]" << std::endl);
                FIXOUTVEC(config::Info,"Initial spin position:",ucm.GetInitS(i,0),ucm.GetInitS(i,1),ucm.GetInitS(i,2));
                FIXOUT(config::Info,"K_1^u: " << ucm.GetK1U(i) << " [Joules]" << std::endl);
                FIXOUTVEC(config::Info,"Direction of anisotropy axis:",ucm.GetK1UDir(i,0),ucm.GetK1UDir(i,1),ucm.GetK1UDir(i,2));
                config::printline(config::Info);
            }
            if(ucm.NumAtomsUnitCell() > 5)
            {
                FIXOUT(config::Log,"Unit cell atom:" << i << std::endl);
                FIXOUT(config::Log,"Element:" << ucm.GetElement(i) << std::endl);
                FIXOUTVEC(config::Log,"Positions [real space]:",ucm.GetCoord(i,0)/static_cast<double>(Nk[0]),ucm.GetCoord(i,1)/static_cast<double>(Nk[1]),ucm.GetCoord(i,2)/static_cast<double>(Nk[2]));
                FIXOUTVEC(config::Log,"Positions [k-lattice]:",ucm.GetCoord(i,0),ucm.GetCoord(i,1),ucm.GetCoord(i,2));
                FIXOUT(config::Log,"Part of sublattice:" << ucm.GetSublattice(i) << std::endl);
                FIXOUT(config::Log,"Damping:" << ucm.GetDamping(i) << std::endl);
                FIXOUT(config::Log,"Gyromagnetic ratio:" << ucm.GetGamma(i) << " [gamma_free]" << std::endl);
                FIXOUT(config::Log,"Magnetic moment:" << ucm.GetMu(i)  << " [muB]" << std::endl);
                FIXOUTVEC(config::Log,"Initial spin position:",ucm.GetInitS(i,0),ucm.GetInitS(i,1),ucm.GetInitS(i,2));
                FIXOUT(config::Log,"K_1^u: " << ucm.GetK1U(i) << " [Joules]" << std::endl);
                FIXOUTVEC(config::Log,"Direction of anisotropy axis:",ucm.GetK1UDir(i,0),ucm.GetK1UDir(i,1),ucm.GetK1UDir(i,2));
                config::printline(config::Log);
            }

        }
        if(ucm.NumAtomsUnitCell()>5)
        {
            for(unsigned int i = 0 ; i < 4 ; i++)
            {
                FIXOUT(config::Info,"             . . . " << "    . . ." << std::endl);
            }
            FIXOUT(config::Info,"FOR COMPLETE UNIT CELL INFORMATION SEE LOG FILE:" << "   log.dat" << std::endl);
            for(unsigned int i = 0 ; i < 4 ; i++)
            {
                FIXOUT(config::Info,"             . . . " << "    . . ." << std::endl);
            }
            config::printline(config::Info);

        }
        //for the use of the interaction matrix we require that each magnetic species
        //(in the unit cell) has the same magnetic moment
        unsigned int errStatus=ucm.CheckSpecies();
        if(errStatus!=0)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage(ucm.errorCode(errStatus));
        }

        //the 5 entries for each spin correspond to
        // 0,1,2 - the x,y,z positions in the unit cell
        // 3 is the magnetic species number
        // 4 is te atom in the unit cell
        lu.resize(nspins,5);
        lu.IFill(0);
        //the 3 bits of information:
        // 0 - the magnetic atom number
        // 1 - the magnetic species type
        coords.resize(zpdim[0]*Nk[0],zpdim[1]*Nk[1],zpdim[2]*Nk[2],2);
        //IF element 0 has the following numbers
        //-2 THEN here corresponds to empty k-mesh point
        //-1 THEN corresponds to an empty k-mesh point but with an imaginary atom


        //there for the determination of the interaction matrix
        coords.IFill(-2);
        std::ofstream sloc("structure.xyz");
        if(sloc.is_open()!=true)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing atomic position information");
        }
        unsigned int atom_counter=0;
        //for counting the number of species
        unsigned int spec_counter[ucm.GetNMS()];
        for(unsigned int i = 0 ; i < ucm.GetNMS() ; i++)
        {
            spec_counter[i]=0;
        }
        //loop over dimensions in x,y and z
        for(unsigned int i = 0 ; i < dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < dim[2] ; k++)
                {
                    for(unsigned int t = 0 ; t < ucm.NumAtomsUnitCell() ; t++)
                    {
                        double ri[3]={0,0,0};
                        coords(i*Nk[0]+ucm.GetCoord(t,0),j*Nk[1]+ucm.GetCoord(t,1),k*Nk[2]+ucm.GetCoord(t,2),0)=atom_counter;
                        coords(i*Nk[0]+ucm.GetCoord(t,0),j*Nk[1]+ucm.GetCoord(t,1),k*Nk[2]+ucm.GetCoord(t,2),1)=ucm.GetSublattice(t);
                        lu(atom_counter,0)=i*Nk[0]+ucm.GetCoord(t,0);
                        lu(atom_counter,1)=j*Nk[1]+ucm.GetCoord(t,1);
                        lu(atom_counter,2)=k*Nk[2]+ucm.GetCoord(t,2);
                        lu(atom_counter,3)=ucm.GetSublattice(t);
                        spec_counter[ucm.GetSublattice(t)]++;
                        lu(atom_counter,4)=t;
                        //1D arrays
                        sublattice[atom_counter]=ucm.GetSublattice(t);
                        gamma[atom_counter]=ucm.GetGamma(t);
                        lambda[atom_counter]=ucm.GetDamping(t);
                        mu[atom_counter]=ucm.GetMu(t);
                        llgpf[atom_counter]=-gamma[atom_counter]/((1.0+lambda[atom_counter]*lambda[atom_counter]));
                        anis::k1u[atom_counter]=2.0*ucm.GetK1U(t)/(llg::muB*mu[atom_counter]);
                        anis::k1udir(atom_counter,0)=ucm.GetK1UDir(t,0);
                        anis::k1udir(atom_counter,1)=ucm.GetK1UDir(t,1);
                        anis::k1udir(atom_counter,2)=ucm.GetK1UDir(t,2);
                        // The real space position is equal to
                        // r_x = ((k_x+i*N_{k,x})/N_{k,x})*lattice constant_x
                        // where k_x is the location of the atom in the unit cell
                        // N_{k,x} is the number of k-points in the x direction
                        // This can of course be generalised the each direction

                        rx[atom_counter]=((ucm.GetCoord(t,0)+static_cast<double>(i*Nk[0]))/static_cast<double>(Nk[0]))*abc[0];
                        ry[atom_counter]=((ucm.GetCoord(t,1)+static_cast<double>(j*Nk[1]))/static_cast<double>(Nk[1]))*abc[1];
                        rz[atom_counter]=((ucm.GetCoord(t,2)+static_cast<double>(k*Nk[2]))/static_cast<double>(Nk[2]))*abc[2];
                        atom_counter++;

                    }
                }
            }
        }
        //a bit of error checking for consistency with number of spins
        if(atom_counter!=nspins)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Number of atoms placed on mesh is more or less than expected.");
        }
        for(unsigned int i = 0 ; i < ucm.GetNMS() ; i++)
        {
            ucm.SetNES(i,spec_counter[i]);
            std::stringstream sstr;
            sstr << "Number of species " << i << ":";
            std::string str=sstr.str();
            FIXOUT(config::Info,str.c_str() << ucm.GetNES(i) << " [ " << 100*(static_cast<double>(ucm.GetNES(i))/static_cast<double>(nspins)) << "% ]" << std::endl);
        }

        if(nms > 5)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Too many magnetic types. If you want more edit cuint.cu variable MAXNSPEC (a #define).\n This is because the constant memory allocation for the species dependent variables cannot be\ndynamic.");
        }
        //sloc << "#This file contains the positions of the magnetic species and their type" << std::endl;
        //sloc << "# Element - x [A] - y [A] - z[A]" << std::endl;
        sloc << nspins << "\n\n";
        for(unsigned int i = 0 ; i < nspins; i++)
        {
            sloc << ucm.GetElement(lu(i,4)) << "\t" << rx[i]/1e-10 << "\t" << ry[i]/1e-10 << "\t" << rz[i]/1e-10 << std::endl;
        }
        sloc.close();
        if(sloc.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close structure files.");
        }
    }
}
