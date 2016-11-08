// File: util_exch_striction.cpp
// Author:Tom Ostler
// Created: 18 Dec 2016
// Last-modified: 19 Oct 2016 14:33:06
// Contains useful functions and classes
// that pertain to exchange striction
#include "../inc/util.h"
#include "../inc/llg.h"
#include "../inc/arrays.h"
#include "../inc/fields.h"
#include "../inc/defines.h"
#include "../inc/spins.h"
#include "../inc/intmat.h"
#include "../inc/geom.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include <string>
#include <sstream>
namespace util
{
    //exch stric no. pairs
    unsigned int esnp,maxuc;
    Array2D<unsigned int> espairs;
    Array3D<int> esuc;
    Array<unsigned int> esnum;
    Array<double> lambda,esP;
    std::ofstream outESP;
    bool escheck;
    void calc_es()
    {
        //zero the esP array that stores the polarisation terms (number of pairs)
        esP.IFill(0);
        //loop over the unit cells
        for(int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(int k = 0 ; k < geom::dim[2] ; k++)
                {
//                    std::cout << i << "\t" << j << "\t" << k << std::endl;
                    //loop over the number of pairs
                    for(unsigned int np = 0 ; np < esnp ; np++)
                    {
                        //get unit cell number of first member of pair
                        unsigned int p1aiuc=espairs(np,0);
                        //get the unit cell number of the second member of pair
                        unsigned int p2aiuc=espairs(np,1);
                        //get the spin id of the first unit cell member of the pair belonging to i,j,k
                        unsigned int p1id=geom::atnolu(i,j,k,p1aiuc);
                        //loop over the number unit cells of the pairs
                        for(unsigned int uc = 0 ; uc < esnum(np) ; uc++)
                        {
                            //get the unit cell lookup
                            int luvec[3]={esuc(np,uc,0),esuc(np,uc,1),esuc(np,uc,2)};
                            int uccoord[3]={i+luvec[0],j+luvec[1],k+luvec[2]};
                            for(unsigned int d = 0 ; d < 3 ; d++)
                            {
                                if(config::pbc[d] && uccoord[d]<0)//then move to the other size
                                {
                                    uccoord[d]=geom::dim[d]+uccoord[d];
                                }
                                else if(config::pbc[d] && uccoord[d]>=geom::dim[d])
                                {
                                    uccoord[d]=uccoord[d]-geom::dim[d];
                                }
                            }
                            //now that we know which unit cell and atom in that unit cell to find
                            //we need to get its spin id
                            unsigned int p2id=geom::atnolu(i,j,k,p2aiuc);
                            //find their spin values
                            double p1s[3]={spins::Sx(p1id),spins::Sy(p1id),spins::Sz(p1id)};
                            double p2s[3]={spins::Sx(p2id),spins::Sy(p2id),spins::Sz(p2id)};
                            //find the dot product
                            double dotp=p1s[0]*p2s[0]+p1s[1]*p2s[1]+p1s[2]*p2s[2];
                            //add to the polarisation array
                            esP(np)+=dotp;
                        }
                        //get the unit cell lookup direction and check if we need to use PBCs
                    }//end of np loop
                }//end of k loop
            }//end of j loop
        }//end of i loop

        //normalise
        for(unsigned int np = 0 ; np < esnp ; np++)
        {
            esP(np)*=(lambda(np)/(static_cast<double>(geom::dim[0]*geom::dim[1]*geom::dim[2]*esnum(np))));
            outESP << "\t" << esP(np);
        }
        outESP << std::endl;


    }
}
