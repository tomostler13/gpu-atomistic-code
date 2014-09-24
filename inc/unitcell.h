// File: array.h
// Author: Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 24 Sep 2014 11:33:20
#ifndef __UNITCELL_H__
#define __UNITCELL_H__
#include "../inc/arrays.h"
#include "../inc/error.h"
#include <vector>
#include <cmath>
//This class contains all of the information about the unit cell
//It does however get converted to more streamlined arrays for high performance
class unitCellMembers
{
    public:
        //default constructor
        //nms - number of magnetic species
        //size - number of atoms in the unit cell
        //elements - string holding element type
        //sublattice - which sublattice (magnetic species type) does this unit cell atom belog to?
        //base_mom - the number of unique moments there should be
        unitCellMembers(): nms(0), size(0), coords(0,0), elements(0), damping(0), mu(0), gamma(0), sublattice(0), initspin(0,0), lambda(0), base_mom(0){}
        //constructor
        unitCellMembers(unsigned int nauc): coords(nauc,3), elements(nauc), damping(nauc), mu(nauc), gamma(nauc), sublattice(nauc), initspin(nauc,3), size(nauc), base_mom(nms) {}
        //destructor
        ~unitCellMembers(){clean();}
        inline void init(unsigned int nauc,unsigned int num_mag)
        {
            //contains the atom coords
            coords.resize(nauc,4);
            //contains the string for the element
            elements.resize(nauc);
            //number of atoms in the unit cell
            size=nauc;
            //damping of each atom in the unit cell
            damping.resize(nauc);
            //gyromagnetic ratio
            gamma.resize(nauc);
            //which sublattice does this belong to
            sublattice.resize(nauc);
            //the intial spin configuration
            initspin.resize(nauc,3);
            //initialse mu
            mu.resize(nauc);
            //set the number of sublattices (num mag species)
            nms=num_mag;
            //resize the base_mom array
            base_mom.resize(nms);
        }
        inline std::string GetElement(unsigned int t)
        {
            assert(t<size);
            return(elements[t]);
        }
        //first 3 arguements are the positions of atom t in unit cell. Second is the atom in the unit cell
        inline void SetPosVec(double x,double y,double z,unsigned int t)
        {
            coords(t,0)=x;
            coords(t,1)=y;
            coords(t,2)=z;
        }
        //first arguement is the string of the atom type. Second is the atom in the unit cell
        inline void SetElement(std::string str,unsigned int t)
        {
            elements[t]=str;
        }
        inline void SetLambda(double x,unsigned int t)
        {
            lambda[t]=x;
        }
        inline void SetGamma(double x,unsigned int t)
        {
            gamma[t]=x;
        }
        inline void SetMu(double x,unsigned int t)
        {
            mu[t]=x;
        }
        inline void SetDamping(double x,unsigned int t)
        {
            damping[t]=x;
        }
        inline void SetSublattice(unsigned int s,unsigned int t)
        {
            assert(s<nms);//if the sublattice number is greater than the number of sublattices that should be there then there is something wrong
            if(s>nms)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("You have specified a given number of sublattices but one of your atoms in your unit cell has a species number greater than this.");
            }
            sublattice[t]=s;
        }
        inline double GetMuBase(unsigned int s)
        {
            return(base_mom[s]);
        }
        unsigned int CheckSpecies()
        {
            //loop over the atoms in the unit cell and determine the first
            //atom with a given sublattice number
            for(unsigned int s = 0 ; s < nms ; s++)
            {
                bool found_base=false;
                for(unsigned int i = 0 ; i < size ; i++)
                {
                    if(sublattice[i]==s)
                    {
                        base_mom[s]=mu[i];
                        found_base=true;
                    }
                    if(found_base)
                    {
                        break;
                    }
                }
            }
            //now we loop over all atoms once more to ensure that all atoms with a
            //given sublattice value have the same magnetic moment (otherwise they
            //should be given another element identifier
            for(unsigned int i = 0 ; i < size ; i++)
            {
                if(fabs(mu[i]-base_mom[sublattice[i]])>1e-12)
                {
                    return(1);
                }
            }
            //if the above loop exits OK then the unit cell moments should be consistent
            return(0);
        }

        inline void SetInitSpin(double sx,double sy,double sz,unsigned int t)
        {
            initspin(t,0)=sx;
            initspin(t,1)=sy;
            initspin(t,2)=sz;
        }
        inline double GetCoord(unsigned int t,unsigned int c)
        {
            return(coords(t,c));
        }
        inline double GetDamping(unsigned int t)
        {
            return(damping(t));
        }
        inline double GetGamma(unsigned int t)
        {
            return(gamma(t));
        }
        inline double GetMu(unsigned int t)
        {
            return(mu(t));
        }
        inline double GetSublattice(unsigned int t)
        {
            return(sublattice(t));
        }
        //first arguement here is the atom in the unit cell. Second is the coordinate (x,y,z)
        inline double GetInitS(unsigned int t,unsigned int c)
        {
            return(initspin(t,c));
        }
        inline void clean()
        {
            coords.clear();
            initspin.clear();
            sublattice.clear();
            mu.clear();
            gamma.clear();
            damping.clear();
        }
        inline unsigned int GetNMS()
        {
            return(nms);
        }
        inline unsigned int NumAtomsUnitCell()
        {
            return(size);
        }
        std::string errorCode(unsigned int errNo)
        {
            if(errNo==0)
            {
                return("No error");
            }
            else if(errNo==1)
            {
                return("A magnetic moment was associated with an atom in the unit cell whose magnetic moment is not consistent with others in that species. Check your unit cell file.");
            }
            else
            {
                return("Error status not recognised");
            }
        }
    private:
        unsigned int size,nms;
        Array2D<double> coords,initspin;
        Array<double> damping,mu,gamma,sublattice,lambda,base_mom;
        std::vector<std::string> elements;
};
#endif /*_UNITCELL_H_*/
