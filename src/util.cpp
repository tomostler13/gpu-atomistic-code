// File: util.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 25 Apr 2013 12:15:59
// Contains useful functions and classes
#include "../inc/util.h"
#include "../inc/arrays.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/intmat.h"
#include "../inc/geom.h"
#include "../inc/error.h"
#include <string>
#include <sstream>
#include <omp.h>
extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}
namespace util
{
    void inverse(double* A, int N)
    {
        int *IPIV = new int[N+1];
        int LWORK = N*N;
        double *WORK = new double[LWORK];
        int INFO;

        dgetrf_(&N,&N,A,&N,IPIV,&INFO);
        dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

        delete [] IPIV;
        delete [] WORK;
    }
    void cpuConvFourier()
    {
        fields::Hkx.IFill(0);
        fields::Hky.IFill(0);
        fields::Hkz.IFill(0);

        //perform convolution in fourier space
        unsigned int i = 0,j = 0, k = 0;
        #pragma omp parallel for private (i)
        for(i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            #pragma omp parallel for private (j)
            for(j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                #pragma omp parallel for private (k) shared (intmat::Hkx,intmat::Hky,intmat::Hkz)
                for(k = 0 ; k < geom::cplxdim ; k++)
                {

                    fields::Hkx(i,j,k)[0]=intmat::Nxx(i,j,k)[0]*spins::Skx(i,j,k)[0]-intmat::Nxx(i,j,k)[1]*spins::Skx(i,j,k)[1]
                        +intmat::Nxy(i,j,k)[0]*spins::Sky(i,j,k)[0]-intmat::Nxy(i,j,k)[1]*spins::Sky(i,j,k)[1]
                        +intmat::Nxz(i,j,k)[0]*spins::Skz(i,j,k)[0]-intmat::Nxz(i,j,k)[1]*spins::Skz(i,j,k)[1];
                    fields::Hkx(i,j,k)[1]=intmat::Nxx(i,j,k)[0]*spins::Skx(i,j,k)[1]+intmat::Nxx(i,j,k)[1]*spins::Skx(i,j,k)[0]
                        +intmat::Nxy(i,j,k)[0]*spins::Sky(i,j,k)[1]+intmat::Nxy(i,j,k)[1]*spins::Sky(i,j,k)[0]
                        +intmat::Nxz(i,j,k)[0]*spins::Skz(i,j,k)[1]+intmat::Nxz(i,j,k)[1]*spins::Skz(i,j,k)[0];

                    fields::Hky(i,j,k)[0]=intmat::Nyx(i,j,k)[0]*spins::Skx(i,j,k)[0]-intmat::Nyx(i,j,k)[1]*spins::Skx(i,j,k)[1]
                        +intmat::Nyy(i,j,k)[0]*spins::Sky(i,j,k)[0]-intmat::Nyy(i,j,k)[1]*spins::Sky(i,j,k)[1]
                        +intmat::Nyz(i,j,k)[0]*spins::Skz(i,j,k)[0]-intmat::Nyz(i,j,k)[1]*spins::Skz(i,j,k)[1];
                    fields::Hky(i,j,k)[1]=intmat::Nyx(i,j,k)[0]*spins::Skx(i,j,k)[1]+intmat::Nyx(i,j,k)[1]*spins::Skx(i,j,k)[0]
                        +intmat::Nyy(i,j,k)[0]*spins::Sky(i,j,k)[1]+intmat::Nyy(i,j,k)[1]*spins::Sky(i,j,k)[0]
                        +intmat::Nyz(i,j,k)[0]*spins::Skz(i,j,k)[1]+intmat::Nyz(i,j,k)[1]*spins::Skz(i,j,k)[0];

                    fields::Hkz(i,j,k)[0]=intmat::Nzx(i,j,k)[0]*spins::Skx(i,j,k)[0]-intmat::Nzx(i,j,k)[1]*spins::Skx(i,j,k)[1]
                        +intmat::Nzy(i,j,k)[0]*spins::Sky(i,j,k)[0]-intmat::Nzy(i,j,k)[1]*spins::Sky(i,j,k)[1]
                        +intmat::Nzz(i,j,k)[0]*spins::Skz(i,j,k)[0]-intmat::Nzz(i,j,k)[1]*spins::Skz(i,j,k)[1];
                    fields::Hkz(i,j,k)[1]=intmat::Nzx(i,j,k)[0]*spins::Skx(i,j,k)[1]+intmat::Nzx(i,j,k)[1]*spins::Skx(i,j,k)[0]
                        +intmat::Nzy(i,j,k)[0]*spins::Sky(i,j,k)[1]+intmat::Nzy(i,j,k)[1]*spins::Sky(i,j,k)[0]
                        +intmat::Nzz(i,j,k)[0]*spins::Skz(i,j,k)[1]+intmat::Nzz(i,j,k)[1]*spins::Skz(i,j,k)[0];

                }
            }
        }
    }
    void copy3vecto1(int size1,double *ia1,double *ia2,double *ia3,double *oa)
    {
        for(int i = 0 ; i < size1 ; i++)
        {
            oa[3*i]=ia1[i];
            oa[3*i+1]=ia2[i];
            oa[3*i+2]=ia3[i];
        }
    }
    void copy3vecto1(int size1,float *ia1,float *ia2,float *ia3,float *oa)
    {
        for(int i = 0 ; i < size1 ; i++)
        {
            oa[3*i]=ia1[i];
            oa[3*i+1]=ia2[i];
            oa[3*i+2]=ia3[i];
        }
    }
    void copy3vecto1(int size1,Array<double> ia1,Array<double> ia2,Array<double> ia3,double *oa)
    {
        for(int i = 0 ; i < size1 ; i++)
        {
            oa[3*i]=ia1[i];
            oa[3*i+1]=ia2[i];
            oa[3*i+2]=ia3[i];
        }
    }
    void copy3vecto1(int size1,Array<float> ia1,Array<float> ia2,Array<float> ia3,float *oa)
    {
        for(int i = 0 ; i < size1 ; i++)
        {
            oa[3*i]=ia1[i];
            oa[3*i+1]=ia2[i];
            oa[3*i+2]=ia3[i];
        }
    }
    void fillfloat(int size,double* da,float* fa)
    {
        for(int i = 0 ; i < size ; i++)
        {
            fa[i]=float(da[i]);
        }
    }
    void fillfloat(int size0,int size1,int size2,Array3D<fftw_complex> da,Array3D<fftwf_complex> fa)
    {
        for(unsigned int i = 0 ; i < size0 ; i++)
        {
            for(unsigned int j = 0 ; j < size1 ; j++)
            {
                for(unsigned int k = 0 ; k < size2 ; k++)
                {
                    fa(i,j,k)[0]=float(da(i,j,k)[0]);
                    fa(i,j,k)[1]=float(da(i,j,k)[1]);
                }
            }
        }
    }
	std::string exec(char* cmd) {
		FILE* pipe = popen(cmd, "r");
		if (!pipe) return "ERROR";
		char buffer[128];
		std::string result = "";
		while(!feof(pipe)) {
			if(fgets(buffer, 128, pipe) != NULL)
				result += buffer;
		}
		pclose(pipe);
		return result;
	}
    void outputSpinsVTU(unsigned int t)
    {
        std::stringstream pvss;
        std::string pv="spinmap";
        pvss << pv << "_" << t << ".vtu";
        std::string pvs = pvss.str();
        std::ofstream pvf(pvs.c_str());
        if(!pvf.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing paraview files");
        }

        pvf << "<?xml version=\"1.0\"?>" << "\n";
        pvf << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
        pvf << "<UnstructuredGrid>" << "\n";
        pvf << "<Piece NumberOfPoints=\""<<geom::nspins<<"\"  NumberOfCells=\"1\">" << "\n";
        pvf << "<PointData Scalar=\"Spin\">" << "\n";
        pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
                    pvf << spins::Sx(i) << "\t" << spins::Sy(i) << "\t" << spins::Sz(i) << "\n";
        }
        pvf << "</DataArray>" << "\n";
        pvf << "</PointData>" << "\n";
        pvf << "<CellData>" << "\n";
        pvf << "</CellData>" << "\n";
        pvf << "<Points>" << "\n";
        pvf << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            pvf << double(geom::lu(i,0))*geom::abc[0] << "\t" << double(geom::lu(i,1))*geom::abc[1] << "\t" << double(geom::lu(i,2))*geom::abc[2] << "\n";
        }

        pvf << "</DataArray>" << "\n";
        pvf << "</Points>" << "\n";
        pvf << "<Cells>" << "\n";
        pvf << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n";
        pvf << "1" << "\n";
        pvf << "</DataArray>" << "\n";
        pvf << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n";
        pvf << "1" << "\n";
        pvf << "</DataArray>" << "\n";
        pvf << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << "\n";
        pvf << "1" << "\n";
        pvf << "</DataArray>" << "\n";
        pvf << "</Cells>" << "\n";
        pvf << "</Piece>" << "\n";
        pvf << "</UnstructuredGrid>" << "\n";
        pvf << "</VTKFile>" << "\n";
        pvf.close();
    }


}

