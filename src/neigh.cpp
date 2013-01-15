// File: neigh.cpp
// Author:Tom Ostler
// Last-modified: 08 Jan 2013 18:41:42
#include "../inc/neigh.h"
#include "../inc/array.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace neigh
{
    bool pbc[3]={false,false,false};
    bool ni;
    //these are the CSR format neighbour lists
    //xadj[i] gives the starting neighbour for neighbour i
    //and ends at (but does not include) xadj[i+1]
    Array<unsigned int> xadj;
    //adjncy[xadj[i]] returns first neighbours
    //.
    //.
    //.
    //adjncy[xadj[i+1]-1] returns last neighbour
    Array<unsigned int> adjncy;
    Array<double> surfArea;

    int bc(int x,int dimx)
    {
        return((x+dimx)%dimx);
    }

    void initNeigh(int argc,char *argv[])
    {
        xadj.resize(geom::ss+1);
        xadj.IFill(0);
        //this will be copied to adjncy but a priori we do not
        //know its size
        std::vector<unsigned int> temp(1);
        std::vector<double> tempsurf(1);
        std::fill(temp.begin(),temp.end(),0);
        std::fill(tempsurf.begin(),tempsurf.end(),0);
        unsigned int nlc=0;
        unsigned int ac=0;
        FIXOUT(config::Info,"Finding neighbours" << std::flush);
        for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                {
                    //if our macrospin exists then find its neighbours
                    if(geom::coords(i,j,k)>-1)
                    {
                        xadj[ac]=nlc;
                        //left neighbour
                        int ln[3]={i-1,j,k};
                        //right neighbour
                        int rn[3]={i+1,j,k};
                        //up neighbour
                        int un[3]={i,j+1,k};
                        //down neighbour
                        int dn[3]={i,j-1,k};
                        //forward neighbour
                        int fn[3]={i,j,k+1};
                        //backward neighbour
                        int bn[3]={i,j,k-1};
                        //boundary conditions for left neighbour
                        if(ln[0]>=0 && ln[0] < geom::dim[0])
                        {
                            if(geom::coords(ln[0],ln[1],ln[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(ln[0],ln[1],ln[2]);
                                //assert(geom::coords(ln[0],ln[1],ln[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[0]*geom::gs[0];
                                nlc++;
                            }
                        }
                        else if(ln[0]<0 && pbc[0]==true)
                        {
                            ln[0]=bc(ln[0],geom::dim[0]);
                            if(geom::coords(ln[0],ln[1],ln[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(ln[0],ln[1],ln[2]);
                            //    assert(geom::coords(ln[0],ln[1],ln[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[0]*geom::gs[0];
                                nlc++;
                            }
                        }
                        //boundary conditions for right neighbour
                        if(rn[0]<geom::dim[0] && rn[0]>=0)
                        {
                            if(geom::coords(rn[0],rn[1],rn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(rn[0],rn[1],rn[2]);
                            //    assert(geom::coords(rn[0],rn[1],rn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[0]*geom::gs[0];
                                nlc++;
                            }
                        }
                        else if(rn[0]>=geom::dim[0] && pbc[0]==true)
                        {
                            rn[0]=bc(rn[0],geom::dim[0]);
                            if(geom::coords(rn[0],rn[1],rn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(rn[0],rn[1],rn[2]);
                            //    assert(geom::coords(rn[0],rn[1],rn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[0]*geom::gs[0];
                                nlc++;
                            }
                        }
                        //boundary conditions for up neighbours
                        if(un[1]<geom::dim[1] && un[1]>=0)
                        {
                            if(geom::coords(un[0],un[1],un[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(un[0],un[1],un[2]);
                            //    assert(geom::coords(un[0],un[1],un[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[1]*geom::gs[1];
                                nlc++;
                            }
                        }
                        else if(un[1]>=geom::dim[1] && pbc[1]==true)
                        {
                            un[1]=bc(un[1],geom::dim[1]);
                            if(geom::coords(un[0],un[1],un[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(un[0],un[1],un[2]);
                            //    assert(geom::coords(un[0],un[1],un[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[1]*geom::gs[1];
                                nlc++;
                            }
                        }
                        //down
                        if(dn[1]>=0 && dn[1]<geom::dim[1])
                        {
                            if(geom::coords(dn[0],dn[1],dn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(dn[0],dn[1],dn[2]);
                            //    assert(geom::coords(dn[0],dn[1],dn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[1]*geom::gs[1];
                                nlc++;
                            }
                        }
                        else if(dn[1]<0 && pbc[1]==true)
                        {
                            dn[1]=bc(dn[1],geom::dim[1]);
                            if(geom::coords(dn[0],dn[1],dn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(dn[0],dn[1],dn[2]);
                            //    assert(geom::coords(dn[0],dn[1],dn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[1]*geom::gs[1];
                                nlc++;
                            }
                        }
                        //forward
                        if(fn[2]>=0 && fn[2] < geom::dim[2])
                        {
                            if(geom::coords(fn[0],fn[1],fn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(fn[0],fn[1],fn[2]);
                            //    assert(geom::coords(fn[0],fn[1],fn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[2]*geom::gs[2];
                                nlc++;
                            }
                        }
                        else if(fn[2]>=geom::dim[2] && pbc[2]==true)
                        {
                            fn[2]=bc(fn[2],geom::dim[2]);
                            if(geom::coords(fn[0],fn[1],fn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(fn[0],fn[1],fn[2]);
                            //    assert(geom::coords(fn[0],fn[1],fn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[2]*geom::gs[2];
                                nlc++;
                            }
                        }
                        //backwards
                        if(bn[2]<geom::dim[2] && bn[2]>=0)
                        {
                            if(geom::coords(bn[0],bn[1],bn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(bn[0],bn[1],bn[2]);
                            //    assert(geom::coords(bn[0],bn[1],bn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[2]*geom::gs[2];
                                nlc++;
                            }
                        }
                        else if(bn[2]<0 && pbc[2]==true)
                        {
                            bn[2]=bc(bn[2],geom::dim[2]);
                            if(geom::coords(bn[0],bn[1],bn[2])>-1)
                            {
                                temp.push_back(1);
                                tempsurf.push_back(1);
                                temp[nlc]=geom::coords(bn[0],bn[1],bn[2]);
                            //    assert(geom::coords(bn[0],bn[1],bn[2])<geom::ss);
                                tempsurf[nlc]=geom::gs[2]*geom::gs[2];
                                nlc++;
                            }
                        }
                        ac++;
                    }
                }
            }
        }
        config::Info << "Done" << std::endl;
        xadj[geom::ss]=nlc;
        FIXOUT(config::Info,"Neighbour list size:" << nlc << std::endl);
        adjncy.resize(nlc);
        surfArea.resize(nlc);
        for(unsigned int i = 0 ; i < nlc ; i++)
        {
            adjncy[i]=temp[i];
            surfArea[i]=tempsurf[i];
        }
        temp.clear();
        tempsurf.clear();

//===========================FOR DEBUGGING=====================================================
        //this can be uncommented and cout redirected to a file (e.g. file.txt)
        //This can be plotted with the numbers as points using
        //gnuplot> plot "file.txt" u 2:3:1 with labels
        /*
        for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                {
                    unsigned int atom=geom::coords(i,j,k);
                    std::stringstream sstr;
                    sstr << "\"" << atom << "\"\t" << i << "\t" << j << "\t" << k << std::endl;
                    std::cout << sstr.str() << std::endl;
                }
            }
        }
//===========================FOR DEBUGGING=====================================================
        */
        //this prints out which cell numbers interact with the others
        /*
        int counter=0;
        for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                {
                    unsigned int at=geom::coords(i,j,k);
                    std::cout << "Cell " << at << " interacts with" << std::endl;
                    for(unsigned int j = xadj[counter] ; j < xadj[counter+1] ; j++)
                    {
                        std::cout << adjncy[j] << "\t" << std::flush;
                    }
                    counter++;
                    std::cout << std::endl;
                }
            }
        }
        */
//=============================================================================================
        //the neighbour list has been initialized
        ni=true;
    }
}
