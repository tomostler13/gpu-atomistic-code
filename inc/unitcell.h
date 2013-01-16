// File: array.h
// Author: Tom Ostler
// Created: 16 Jan 2013
// Last-modified: 16 Jan 2013 17:46:49
#ifndef __UNITCELL_H__
#define __UNITCELL_H__
#include "../inc/arrays.h"
class unitCellMembers
{
    public:
        //default constructork
        unitCellMembers(): data(0,0) {}
        //constructor
        unitCellMembers(unsigned int nauc): data(nauc,4) {}
        //destructor
        ~unitCellMembers(){clean();}
        inline void init(unsigned int nauc)
        {
            data.resize(nauc,4);
        }
        inline void SetPosVec(double x,double y,double z,unsigned int t)
        {
            data(t,0)=x;
            data(t,1)=y;
            data(t,2)=z;
            data(t,3)=t;
        }
        inline void SetX(double x,unsigned int t)
        {
            data(t,0)=x;
        }
        inline void SetY(double y,unsigned int t)
        {
            data(t,1)=y;
        }
        inline void SetZ(double z,unsigned int t)
        {
            data(t,2)=z;
        }
        inline void SetT(unsigned int t)
        {
            data(t,3)=t;
        }
        inline double GetX(unsigned int t)
        {
            return(data(t,0));
        }
        inline double GetY(unsigned int t)
        {
            return(data(t,1));
        }
        inline double GetZ(unsigned int t)
        {
            return(data(t,2));
        }
        inline double GetComp(unsigned int t,unsigned int comp)
        {
            return(data(t,comp));
        }
        inline void clean()
        {
            data.clear();
        }
    private:
        unsigned int size;
        Array2D<double> data;
};
#endif /*_UNITCELL_H_*/
