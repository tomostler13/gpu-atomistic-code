// File: mat.cpp
// Author:Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 22 Apr 2013 14:41:51
#include "../inc/mat.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/arrays.h"
#include "../inc/random.h"
#include <libconfig.h++>
#include <cassert>
#include <sstream>
#include <string>
#include <cmath>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace mat
{
    //initially set all spins to same damping
    double ilambda=0.0;
    //bath coupling (no units)
    Array<double> lambda;
    //Gyromagnetic ratio
    double gamma=0.0;
    //Bohr magneton
    double muB=9.27e-24;
    //number of species
    unsigned int nspec=1;
    //magnetic moment (muB)
    Array<double> mu,mustore;
    //Concentration
    Array<double> conc;
    //Holds species number
    Array<unsigned int> speclist;
    //thermal term prefactor
    Array<double> sigma;
    std::string place;
    void initMat(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Material details***" << std::endl;
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

        libconfig::Setting &setting = config::cfg.lookup("mat");

        setting.lookupValue("lambda",ilambda);
        FIXOUT(config::Info,"Coupling constant (lambda):" << ilambda << std::endl);
        setting.lookupValue("gamma",gamma);
        FIXOUT(config::Info,"Gyromagnetic ratio:" << gamma << " T^-1 s^-1\n");
        muB=9.27e-24;
        setting.lookupValue("NumSpec",nspec);
        if(config::useintmat==true && nspec > 1)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You cannot have more than 1 species if you want to use the interaction matrix");
        }
        FIXOUT(config::Info,"Number of species:" << nspec << std::endl);
        mustore.resize(nspec);
        for(unsigned int i = 0 ; i < nspec ; i++)
        {
            mustore[i]=setting["mu"][i];
            FIXOUT(config::Info,"Magnetic moments (species/moment):" << i << "  /  " << mustore[i] << std::endl);
        }
        mu.resize(geom::nspins);
        lambda.resize(geom::nspins);
        speclist.resize(geom::nspins);
        sigma.resize(geom::nspins);
        sigma.IFill(0);
        if(nspec>1)
        {
            setting.lookupValue("Place",place);
            FIXOUT(config::Info,"Placing method:" << place << std::endl);
            if(place=="random")
            {
                double sum=0.0;
                Array<int> naes(nspec);

                for(unsigned int i = 0 ; i < nspec ; i++)
                {
                    conc[i]=setting["conc"][i];
                    sum+=conc[i];
                    FIXOUT(config::Info,"Desired concentration (species/concentration):" << i << "  /  " << conc[i] << std::endl);
                    naes[i]=int(geom::nspins*conc[i]);
                    FIXOUT(config::Info,"Number of atoms:" << naes[i] << std::endl);
                    FIXOUT(config::Info,"Actual concentration:" << double(naes[i])/(double(geom::nspins)));
                }
                if(fabs(1-sum)>1e-10)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Concentration not correct");
                }
                for(unsigned int n = 0 ; n < nspec-1 ; n++)
                {

                    bool check=false;
                    unsigned int counter=0;
                    while(check==false)
                    {
                        for(unsigned int i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
                        {
                            for(unsigned int j = 0 ; j < geom::dim[1]*geom::Nk[1] ; j++)
                            {
                                for(unsigned int k = 0 ; k < geom::dim[2]*geom::Nk[2] ; k++)
                                {
                                    //generate a random number between 0 and 1
                                    const double rn=Random::rand();
                                    //check
                                    //      -if rn is less than desired percentage
                                    //      -if an atom exists there
                                    //      -if we are replacing species 0. We do not want to replace on we have just placed
                                    if(rn < conc[n] && geom::coords(i,j,k,0)>-1 && geom::coords(i,j,k,2)==0)
                                    {
                                        geom::coords(i,j,k,2)=n;
                                        speclist(geom::coords(i,j,k,0))=geom::coords(i,j,k,2);
                                        counter++;
                                        if(counter=naes[n])
                                        {
                                            check=true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
