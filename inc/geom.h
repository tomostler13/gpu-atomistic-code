// File: geom.h
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 16 Jan 2013 15:53:39
#include "../inc/arrays.h"
#include <string>
#ifndef _GEOM_H_
#define _GEOM_H_
namespace geom
{
    extern unsigned int dim[],nauc;
    void initGeom(int argc,char *argv[]);
    extern Array2D<double> L,Linv;
    extern Array<double> abc;
    extern Array<unsigned int> Nk;
    class unitCellMembers
    {
        public:
            //constructor
            unitCellMembers(unsigned int nauc): data(nauc,4) {}
            //destructor
            ~unitCellMembers(){clean();}
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

}
#endif /*_GEOM_H_*/
