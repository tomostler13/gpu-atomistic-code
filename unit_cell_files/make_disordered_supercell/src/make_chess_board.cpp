// File: make_chess_board.cpp
// Author:Tom Ostler
// Created: 22 Nov 2014
<<<<<<< HEAD
// Last-modified: 02 Dec 2014 16:50:06
=======
// Last-modified: 02 Dec 2014 19:43:34
>>>>>>> 325eed9487ec64f9ab051b54d5d885be0c06afed

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
#include "../../../inc/random.h"
int main(int argc,char *argv[])
{
    //for outputting info
    std::ofstream Info("structure_info.dat");
    unsigned int seed=0,nspins=0,nspinsperarea;
    unsigned int dim[3]={0,0,0},globdim[3]={0,0,0},nk[3]={1,1,1},nms = 1,nauc=1;
    Array<double> damp,gamma,anis,mu;
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
    errstatus=setting.lookupValue("seed",seed);
    if(errstatus)
    {
        if(seed>0)
        {
            FIXOUT(Info,"Random number seed:" << seed << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The seed must be greater than 0");
        }
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not read random number generator seed.");
    }
    Random::seed(seed,seed+100);

    errstatus=setting.lookupValue("nms",nms);

    Info << "---------------------------------------------------------------" << std::endl;
    Info << "            Magnetic Species Info                              " << std::endl;
    if(errstatus)
    {
        FIXOUT(Info,"Number of magnetic species:" << nms << std::endl);

        amounts.resize(2,nms);amounts.IFill(0);
        damp.resize(nms);damp.IFill(0);
        mu.resize(nms);mu.IFill(0);
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
        mu[i]=setting["mu"][i];
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
        FIXOUT(Info,"Moment (mu_B):" << mu[i] << std::endl);
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
    unsigned int dimension=0;
    for(unsigned int i = 0 ; i < 3 ; i++)
    {
        if(globdim[i]>1)
        {
            dimension++;
        }
    }
    FIXOUT(Info,"Dimension of \"chess board\":" << dimension << std::endl);
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
    nspins=atom_counter;
    nspinsperarea=static_cast<int>(static_cast<double>(nspins)/static_cast<double>(globdim[0]*globdim[1]*globdim[2])+0.5);
    FIXOUT(Info,"Total number of spins:" << nspins << std::endl);
    FIXOUT(Info,"Number of spins per volume:" << nspinsperarea << std::endl);
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

    std::ofstream compmap("compmap.dat");
    if(!compmap.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open compmap.dat file.");
    }
    compmap << "# Cell in x - Cell in y - Cell in z - Desired \% for species 1 ...... - Actual \% for species 1 ....... " << std::endl;
    compmap << "# note that the \% of species 0 is inferred from the amounts of species 1, 2 ...." << std::endl;
    //loop over the unit cells in each direction
    for(unsigned int i = 0 ; i < globdim[0] ; i++)
    {
        unsigned int xoffs=i*dim[0]*nk[0];
        for(unsigned int j = 0 ; j < globdim[1] ; j++)
        {
            unsigned int yoffs=j*dim[1]*nk[1];
            for(unsigned int k = 0 ; k < globdim[2] ; k++)
            {
                unsigned int zoffs=k*dim[2]*nk[2];
                //lookup which composition we want
                unsigned int comp=((i+j+k)%nms);

                //std::cout << i << "\t" << j << "\t" << k << "\t" << comp << std::endl;
                compmap << i << "\t" << j << "\t" << k;
                //loop over each species after 0 (i.e. the second and third etc) until the desired percentage
                for(unsigned int species = 1 ; species < nms ; species++)
                {
                    double percent=0.0;
                    unsigned int countamount=0;
                    bool placed=false;
                    while(placed==false)
                    {
                        //now loop over the k-points within this "chess square"
                        for(unsigned int ii = 0 ; ii < dim[0]*nk[0] ; ii++)
                        {
                            unsigned int xcoord=xoffs+ii;
                            for(unsigned int jj = 0 ; jj < dim[1]*nk[1] ; jj++)
                            {
                                unsigned int ycoord=yoffs+jj;
                                for(unsigned int kk = 0 ; kk < dim[2]*nk[2] ; kk++)
                                {
                                    unsigned int zcoord=zoffs+kk;

                                    //std::cout << "Coord\t" << xcoord << "\t" << ycoord << "\t" << zcoord << "\t" << coords(xcoord,ycoord,zcoord) << std::endl;
                                    //we only want to replace atoms that should
                                    //exist (>-1) and ones that are not species
                                    //that we have placed after placing species 0
                                    if(coords(xcoord,ycoord,zcoord)==0)
                                    {
                                        if(Random::rand()*100.0 <= amounts(comp,species))
                                        {
                                            //then place an atom
                                            coords(xcoord,ycoord,zcoord)=species;
                                            countamount++;
                                            percent=(static_cast<double>(countamount)/static_cast<double>(nspinsperarea))*100.0;
                                            //std::cout << xcoord << "\t" << ycoord << "\t" << zcoord << "\tpercent = " << percent << std::endl;
                                            if(percent >= amounts(comp,species))
                                            {
                                                placed=true;
                                            }
                                        }
                                    }
                                    if(placed){break;}
                                }
                                if(placed){break;}
                            }
                            if(placed){break;}
                        }
                    }

                    compmap << "\t" << amounts(comp,species) << "\t" << percent;
                }
                compmap << std::endl;

            }
        }
    }
    std::ofstream ucf("disordered.ucf");
    if(!ucf.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open unit cell file.");
    }
    ucf << nms << std::endl;
    ucf << nspins << std::endl;

    const double sizes[3]={static_cast<double>(dim[0]*globdim[0]*nk[0]),
                           static_cast<double>(dim[1]*globdim[1]*nk[1]),
                           static_cast<double>(dim[2]*globdim[2]*nk[2])};
    for(unsigned int species = 0 ; species < nms ; species++)
    {
/*        std::stringstream sstr;
        sstr << "spinmap_spec_" << species << ".dat";
        std::string str=sstr.str();
        std::ofstream spinmap(str.c_str());
        if(!spinmap.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing spinmap.");
        }
*/
        for(unsigned int i = 0 ; i < dim[0]*globdim[0]*nk[0] ; i++)
        {
            for(unsigned int j = 0 ; j < dim[1]*globdim[1]*nk[1] ; j++)
            {
                for(unsigned int k = 0 ; k < dim[2]*globdim[2]*nk[2] ; k++)
                {
                    if(coords(i,j,k)==species)
                    {
                        //spinmap << i << "\t" << j << "\t" << k << std::endl;
                        ucf << species << "\t" << static_cast<double>(i)/sizes[0] << "\t" << static_cast<double>(j)/sizes[1] << "\t" << static_cast<double>(k)/sizes[2] << "\t" << mu[species] << "\t" << damp[species] << "\t" << gamma[species] << "\t" << spec[species] << "\t" << spin(species,0) << "\t" << spin(species,1) << "\t" << spin(species,2) << "\t" << anis[species] << "\t" << anisvec(species,0) << "\t" << anisvec(species,1) << "\t" << anisvec(species,2) << std::endl;

                    }
                }
            }
        }
/*        spinmap.close();
        if(spinmap.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close spinmap file.");
<<<<<<< HEAD
        }
    }*/
=======
        }*/
>>>>>>> 325eed9487ec64f9ab051b54d5d885be0c06afed
    }
    ucf.close();
    if(ucf.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close the unit cell file.");
    }
    compmap.close();
    return(EXIT_SUCCESS);
}
