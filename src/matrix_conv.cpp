//File matrix_conv.cpp
// Author: Tom Ostler
// Created: 07 Oct 2014
// Last-modified: 15 Jun 2016 11:31:33
// The routines within this file convert a 2D matrix to
// a number of formats depending on the routine used. The
// return structures depend on the storage format
#include "../inc/arrays.h"
#include "../inc/matrix_conv.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include <string>
#include <iostream>
#include <cstdlib>
namespace matconv
{
    //First arguement is number of spins
    //The second two arguements here make up the CSR storage
    //diagoffset is an output array that holds the offsets
    //the diagonals that have non-zero elements have been
    //precalculated during the formation of the CSR storage
    //and are stored in checkdiag which we will use to
    //determine the offsets. The numdiag has already been
    //determined. The dataxx, datayy and datazz on input
    //contain the CSR format exchange info. On output they
    //should be the data for the DIA format (they will need
    //be resized)
    void csr_to_dia_diag(unsigned int& N,Array<unsigned int>& xadj,Array<unsigned int>& adjncy,
            Array<int>& offset,unsigned int numdiag,
            Array<int>& checkdiag,
            Array<double>& dataxx,
            Array<double>& datayy,
            Array<double>& datazz
            )
    {
        //resize the offset
        offset.resize(numdiag);
        //temporary arrays to store the dia format data before
        //resizing and copying back to the original arrays
        Array<double> tdataxx,tdatayy,tdatazz;
        tdataxx.resize(numdiag*N);
        tdatayy.resize(numdiag*N);
        tdatazz.resize(numdiag*N);
        //initialise to zero
        for(unsigned int i = 0 ; i < numdiag*N ; i++)
        {
            tdataxx[i]=0;
            tdatayy[i]=0;
            tdatazz[i]=0;
        }
        unsigned int count=0;
        for(unsigned int i = 0 ; i < 2*N-1 ; i++)
        {
            if(checkdiag[i]>0)
            {//then this diagonal has some data in it

                //write the diagonal to the offset
                offset[count]=static_cast<int>(i)-(static_cast<int>(N));
                count++;
            }
        }
        if(count!=numdiag)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The number of diagonals that we have recounted is not correct.");
        }

        //loop over the diagonals
        for(unsigned int i = 0 ; i < numdiag ; i++)
        {
            int diag = offset[i];
            int num_elem=0,num_zeros=0;
            if(diag<=0)
            {
                num_elem=diag+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-diag;
            }
            //so num_elem+num_zeroes = N
            num_zeros=N-num_elem;
            unsigned int istart=0,jstart=0;
            if(diag<=0)
            {
                istart=-diag;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=diag;
            }

            int iloop=istart,jloop=jstart;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {
                //starting value of CSR to lookup
                unsigned int start=xadj[iloop],end=xadj[iloop+1];
                for(unsigned int k = start ; k < end ; k++)
                {
                    //now need to determine if k is on our diagonal
                    if(adjncy[k]==jloop)
                    {
                        //then we have found a neighbour that is on our diagonal
                        if(diag<0)
                        {
                            tdataxx[i*N+num_zeros+j]=dataxx[k];
                            tdatayy[i*N+num_zeros+j]=datayy[k];
                            tdatazz[i*N+num_zeros+j]=datazz[k];
                        }
                        else
                        {
                            tdataxx[i*N+j]=dataxx[k];
                            tdatayy[i*N+j]=datayy[k];
                            tdatazz[i*N+j]=datazz[k];
                        }
                    }
                }

                iloop++;
                jloop++;
            }
        }
        //we no longer want the data in dataxx, datayy and datazz (well
        //we want the DIA data in there).
        dataxx.resize(numdiag*N);
        datayy.resize(numdiag*N);
        datazz.resize(numdiag*N);
        for(unsigned int i = 0 ; i < numdiag*N ; i++)
        {
            dataxx[i]=tdataxx[i];
            datayy[i]=tdatayy[i];
            datazz[i]=tdatazz[i];
        }
        //this is a bit pointless because they are going out of scope
        tdataxx.clear();
        tdatayy.clear();
        tdatazz.clear();

    }
    void csr_to_dia_offdiag(unsigned int& N,Array<unsigned int>& xadj,Array<unsigned int>& adjncy,
            Array<int>& offset,unsigned int numdiag,
            Array<int>& checkdiag,
            Array<double>& dataxy,
            Array<double>& dataxz,
            Array<double>& datayx,
            Array<double>& datayz,
            Array<double>& datazx,
            Array<double>& datazy
            )
    {
        //temporary arrays to store the dia format data before
        //resizing and copying back to the original arrays
        Array<double> tdataxy,tdataxz,tdatayx,tdatayz,tdatazx,tdatazy;
        tdataxy.resize(numdiag*N);
        tdataxz.resize(numdiag*N);
        tdatayx.resize(numdiag*N);
        tdatayz.resize(numdiag*N);
        tdatazx.resize(numdiag*N);
        tdatazy.resize(numdiag*N);
        //initialise to zero
        for(unsigned int i = 0 ; i < numdiag*N ; i++)
        {
            tdataxy[i]=0;
            tdataxz[i]=0;
            tdatayx[i]=0;
            tdatayz[i]=0;
            tdatazx[i]=0;
            tdatazy[i]=0;
        }
        //loop over the diagonals
        for(unsigned int i = 0 ; i < numdiag ; i++)
        {
            int diag = offset[i];
            int num_elem=0,num_zeros=0;
            if(diag<=0)
            {
                num_elem=diag+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-diag;
            }
            //so num_elem+num_zeroes = N
            num_zeros=N-num_elem;
            unsigned int istart=0,jstart=0;
            if(diag<=0)
            {
                istart=-diag;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=diag;
            }

            int iloop=istart,jloop=jstart;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {
                //starting value of CSR to lookup
                unsigned int start=xadj[iloop],end=xadj[iloop+1];
                for(unsigned int k = start ; k < end ; k++)
                {
                    //now need to determine if k is on our diagonal
                    if(adjncy[k]==jloop)
                    {
                        //then we have found a neighbour that is on our diagonal
                        if(diag<0)
                        {
                            tdataxy[i*N+num_zeros+j]=dataxy[k];
                            tdataxz[i*N+num_zeros+j]=dataxz[k];
                            tdatayx[i*N+num_zeros+j]=datayx[k];
                            tdatayz[i*N+num_zeros+j]=datayz[k];
                            tdatazx[i*N+num_zeros+j]=datazx[k];
                            tdatazy[i*N+num_zeros+j]=datazy[k];
                        }
                        else
                        {
                            tdataxy[i*N+j]=dataxy[k];
                            tdataxz[i*N+j]=dataxz[k];
                            tdatayx[i*N+j]=datayx[k];
                            tdatayz[i*N+j]=datayz[k];
                            tdatazx[i*N+j]=datazx[k];
                            tdatazy[i*N+j]=datazy[k];
                        }
                    }
                }

                iloop++;
                jloop++;
            }
        }
        //we no longer want the data in dataxx, datayy and datazz (well
        //we want the DIA data in there). We are going to do these in pairs
        //incase we malloc too much memory
        dataxy.resize(numdiag*N);
        dataxz.resize(numdiag*N);
        for(unsigned int i = 0 ; i < numdiag*N ; i++)
        {
            dataxy[i]=tdataxy[i];
            dataxz[i]=tdataxz[i];
        }
        tdataxy.clear();tdataxz.clear();
        datayx.resize(numdiag*N);
        datayz.resize(numdiag*N);
        for(unsigned int i = 0 ; i < numdiag*N ; i++)
        {
            datayx[i]=tdatayx[i];
            datayz[i]=tdatayz[i];
        }
        tdatayx.clear();
        tdatayz.clear();
        datazx.resize(numdiag*N);
        datazy.resize(numdiag*N);
        for(unsigned int i = 0 ; i < numdiag*N ; i++)
        {
            datazx[i]=tdatazx[i];
            datazy[i]=tdatazy[i];
        }
        tdatazx.clear();tdatazy.clear();
    }

    // This calculates the diagonal offsets for the DIA format
    // and determines the data.
    // The first arguement is the array to store the offsets
    // The second is the J matrix (4 dimensions because it
    // represents all elements of the interactio matrix).
    // This function is overloaded. If the arguement list
    // contains only one offset array then we calculate that
    // for the diagonals of the exchange tensor.
    void conv_intmat_to_dia(Array<int>& offset,Array4D<double>& JMat,unsigned int& num_diag,unsigned int& nN,Array<double>& dataxx,Array<double>& datayy,Array<double>& datazz)
    {
        //for counting the number of non-zero diagonals
        unsigned int dia_count=0;
        int N=static_cast<int>(nN);
        //set the number of diagonals to zero
        num_diag=0;
        for(int i = -N+1 ; i < N ; i++)
        {
            //calculate the number of elements in the diagonal
            int num_elem=0;
            if(i<=0)
            {
                num_elem=i+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-i;
            }
            unsigned int istart=0,jstart=0;
            if(i<=0)
            {
                istart=-i;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=i;
            }
            int iloop=istart,jloop=jstart;

            double abs_J=0.0;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {
                abs_J+=(fabs(JMat(iloop,jloop,0,0))+fabs(JMat(iloop,jloop,1,1))+fabs(JMat(iloop,jloop,2,2)));
                iloop++;
                jloop++;
            }
            if(abs_J>1e-30)
            {
                num_diag++;
            }
        }
        //resize the offset array
        offset.resize(num_diag);
        dataxx.resize(num_diag*N);
        datayy.resize(num_diag*N);
        datazz.resize(num_diag*N);
        dataxx.IFill(0);datayy.IFill(0);datazz.IFill(0);
        offset.IFill(0);
        unsigned int nd=0;
        //loop over the diagonals again and set the value of the offset
        for(int i = -N+1 ; i < N ; i++)
        {
            //calculate the number of elements in the diagonal (that are not zero by default)
            int num_elem=0,num_zeros=0;
            if(i<=0)
            {
                num_elem=i+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-i;
            }
            //so num_elem+num_zeroes = N
            num_zeros=N-num_elem;
            unsigned int istart=0,jstart=0;
            if(i<=0)
            {
                istart=-i;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=i;
            }
            int iloop=istart,jloop=jstart;


            double abs_J=0.0;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {

                abs_J+=(fabs(JMat(iloop,jloop,0,0))+fabs(JMat(iloop,jloop,1,1))+fabs(JMat(iloop,jloop,2,2)));
                iloop++;
                jloop++;
            }
            //if there was a single element in that diagonal that was greater than zero
            if(abs_J>1e-30)
            {
                iloop=istart;
                jloop=jstart;
                for(int j = 0 ; j < num_elem ; j++)
                {
                    if(i<0)
                    {
                        dataxx[N*nd+num_zeros+j]=JMat(iloop,jloop,0,0);
                        datayy[N*nd+num_zeros+j]=JMat(iloop,jloop,1,1);
                        datazz[N*nd+num_zeros+j]=JMat(iloop,jloop,2,2);
                    }
                    else
                    {
                        dataxx[N*nd+j]=JMat(iloop,jloop,0,0);
                        datayy[N*nd+j]=JMat(iloop,jloop,1,1);
                        datazz[N*nd+j]=JMat(iloop,jloop,2,2);
                    }
                    iloop++;
                    jloop++;
                }
                offset[nd]=i;
                nd++;
            }
        }
        //for debugging the data array
/*        for(int i = 0 ; i < num_diag ; i++)
        {
            int os=offset[i];
            std::cout << "Offset " << i << " is " << os << std::endl;
            for(int j = 0 ; j < N ; j++)
            {
                std::cout << j << "\t" << i*num_diag+j << "\t" << datazz[i*num_diag+j] << std::endl;
            }
            std::cin.get();
        }*/

        if(nd!=num_diag)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The number of diagonals upon checking was not consisten");
        }
    }
    void conv_intmat_to_dia(Array<int>& diagoffset,Array<int>& offdiagoffset,Array4D<double>& JMat,unsigned int& num_diag,unsigned int& num_off_diag,unsigned int& nN,Array<double>& dataxx,Array<double>& dataxy,Array<double>& dataxz,Array<double>& datayx,Array<double>& datayy,Array<double>& datayz,Array<double>& datazx,Array<double>& datazy,Array<double>& datazz)
    {
        //for counting the number of non-zero diagonals
        unsigned int dia_count=0;
        int N=static_cast<int>(nN);
        //set the number of diagonals to zero
        num_diag=0,num_off_diag=0;
        for(int i = -N+1 ; i < N ; i++)
        {
            //calculate the number of elements in the diagonal
            int num_elem=0;
            if(i<=0)
            {
                num_elem=i+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-i;
            }
            unsigned int istart=0,jstart=0;
            if(i<=0)
            {
                istart=-i;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=i;
            }
            int iloop=istart,jloop=jstart;

            double abs_J=0.0;
            double abs_offJ=0.0;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {
                abs_J+=fabs(JMat(iloop,jloop,0,0))+fabs(JMat(iloop,jloop,1,1))+fabs(JMat(iloop,jloop,2,2));
                abs_offJ+=fabs(JMat(iloop,jloop,0,1))+fabs(JMat(iloop,jloop,0,2))+fabs(JMat(iloop,jloop,1,0))+fabs(JMat(iloop,jloop,1,2))+fabs(JMat(iloop,jloop,2,0))+fabs(JMat(iloop,jloop,2,1));
                iloop++;
                jloop++;
            }
            if(abs_J>1e-30)
            {
                num_diag++;
            }
            if(abs_offJ>1e-30)
            {
                num_off_diag++;
            }
        }
        //resize the offset array
        diagoffset.resize(num_diag);
        offdiagoffset.resize(num_off_diag);
        dataxx.resize(num_diag*N);
        datayy.resize(num_diag*N);
        datazz.resize(num_diag*N);
        if(num_off_diag>0)
        {
            dataxy.resize(num_off_diag*N);
            dataxz.resize(num_off_diag*N);
            datayx.resize(num_off_diag*N);
            datayz.resize(num_off_diag*N);
            datazx.resize(num_off_diag*N);
            datazy.resize(num_off_diag*N);
            dataxy.IFill(0);dataxz.IFill(0);datayx.IFill(0);datayz.IFill(0);datazx.IFill(0);datazy.IFill(0);
            dataxx.IFill(0);datayy.IFill(0);datazz.IFill(0);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You have specified that you want to include the off diagonal terms in the exchange matrix, however their values are zero.");
        }
        diagoffset.IFill(0);
        offdiagoffset.IFill(0);
        unsigned int nd=0,nod=0;
        //loop over the diagonals again and set the value of the offset
        for(int i = -N+1 ; i < N ; i++)
        {
            //calculate the number of elements in the diagonal
            int num_elem=0,num_zeros=0;
            if(i<=0)
            {
                num_elem=i+static_cast<int>(N);
            }
            else
            {
                num_elem=static_cast<int>(N)-i;
            }
            num_zeros=N-num_elem;
            unsigned int istart=0,jstart=0;
            if(i<=0)
            {
                istart=-i;
                jstart=0;
            }
            else
            {
                istart=0;
                jstart=i;
            }
            int iloop=istart,jloop=jstart;


            double abs_J=0.0;
            double abs_offJ=0.0;
            for(unsigned int j = 0 ; j < num_elem ; j++)
            {

                abs_J+=fabs(JMat(iloop,jloop,0,0))+fabs(JMat(iloop,jloop,1,1))+fabs(JMat(iloop,jloop,2,2));
                abs_offJ+=fabs(JMat(iloop,jloop,0,1))+fabs(JMat(iloop,jloop,0,2))+fabs(JMat(iloop,jloop,1,0))+fabs(JMat(iloop,jloop,1,2))+fabs(JMat(iloop,jloop,2,0))+fabs(JMat(iloop,jloop,2,1));
                iloop++;
                jloop++;
            }
            if(abs_J>1e-30)
            {
                iloop=istart;
                jloop=jstart;
                for(unsigned int j = 0 ; j < num_elem ; j++)
                {
                    if(i<0)
                    {
                        dataxx[num_diag*nd+(num_zeros)+j]=JMat(iloop,jloop,0,0);
                        datayy[num_diag*nd+(num_zeros)+j]=JMat(iloop,jloop,1,1);
                        datazz[num_diag*nd+(num_zeros)+j]=JMat(iloop,jloop,2,2);
                    }
                    else
                    {
                        dataxx[num_diag*nd+j]=JMat(iloop,jloop,0,0);
                        datayy[num_diag*nd+j]=JMat(iloop,jloop,1,1);
                        datazz[num_diag*nd+j]=JMat(iloop,jloop,2,2);
                    }
                    iloop++;
                    jloop++;
                }
                diagoffset[nd]=i;
                nd++;
            }
            if(abs_offJ>1e-30)
            {
                iloop=istart;
                jloop=jstart;
                for(unsigned int j = 0 ; j < num_elem ; j++)
                {
                    if(i<0)
                    {
                        dataxy[N*nod+(num_zeros)+j]=JMat(iloop,jloop,0,1);
                        dataxz[N*nod+(num_zeros)+j]=JMat(iloop,jloop,0,2);
                        datayx[N*nod+(num_zeros)+j]=JMat(iloop,jloop,1,0);
                        datayz[N*nod+(num_zeros)+j]=JMat(iloop,jloop,1,2);
                        datazx[N*nod+(num_zeros)+j]=JMat(iloop,jloop,2,0);
                        datazy[N*nod+(num_zeros)+j]=JMat(iloop,jloop,2,1);
                    }
                    else
                    {
                        dataxy[N*nod+j]=JMat(iloop,jloop,0,1);
                        dataxz[N*nod+j]=JMat(iloop,jloop,0,2);
                        datayx[N*nod+j]=JMat(iloop,jloop,1,0);
                        datayz[N*nod+j]=JMat(iloop,jloop,1,2);
                        datazx[N*nod+j]=JMat(iloop,jloop,2,0);
                        datazy[N*nod+j]=JMat(iloop,jloop,2,1);
                    }
                    iloop++;
                    jloop++;
                }
                offdiagoffset[nod]=i;
                nod++;
            }
        }
        //for debugging the data array
        /*
        for(unsigned int i = 0 ; i < num_diag ; i++)
        {
            int os=offset[i];
            std::cout << "Offset " << i << " is " << os << std::endl;
            for(unsigned int j = 0 ; j < N ; j++)
            {
                std::cout << j << "\t" << data[i*num_diag+j] << std::endl;
            }
            std::cin.get();
        }*/
        if(nd!=num_diag || nod!=num_off_diag)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The number of diagonals upon checking was not consisten");
        }
    }

}
