// File: util.cpp
// Author:Tom Ostler
// Created: 15 Jan 2013
// Last-modified: 21 Mar 2019 01:20:00 PM
// Contains useful functions and classes
#include "../inc/util.h"
#include "../inc/llg.h"
#include "../inc/arrays.h"
#include "../inc/fields.h"
#include "../inc/spins.h"
#include "../inc/intmat.h"
#include "../inc/geom.h"
#include "../inc/error.h"
#include <string>
#include <sstream>
/*extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}*/
namespace util
{
    std::ofstream ofs,sofs;
    Array2D<double> magx,magy,magz,magEachAtUC;
    double lx1,ly1,lz1,lx2,ly2,lz2;
    Array3D<double> mag_species_x,mag_species_y,mag_species_z;
    Array<unsigned int> magDiscSize;
    Array4D<double> magdisc;
    std::string disOutForm;
    unsigned int maxcx=0,maxcy=0,maxcz=0;
    Array3D<unsigned int> nd;
    Array2D<double> nspl;//number of spins per layer
/*    void inverse(double* A, int N)
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
    */
    void cpuConvFourier()
    {
        fields::Hk.IFill(0);

        //perform convolution in fourier space
        register unsigned int i = 0,j = 0, k = 0, s1 = 0, s2=0, alpha=0, beta=0;
        for(i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    for(s1 = 0 ; s1 < geom::ucm.GetNMS() ; s1++)
                    {
                        for(s2 = 0 ; s2 < geom::ucm.GetNMS() ; s2++)
                        {
                            for(alpha = 0 ; alpha < 3 ; alpha++)
                            {
                                for(beta = 0 ; beta < 3 ; beta++)
                                {
                                    fields::Hk(s1,alpha,i,j,k)[0]+=(intmat::Nkab(s1,s2,alpha,beta,i,j,k)[0]*spins::Sk(s2,beta,i,j,k)[0]-intmat::Nkab(s1,s2,alpha,beta,i,j,k)[1]*spins::Sk(s2,beta,i,j,k)[1]);
                                    fields::Hk(s1,alpha,i,j,k)[1]+=(intmat::Nkab(s1,s2,alpha,beta,i,j,k)[0]*spins::Sk(s2,beta,i,j,k)[1]+intmat::Nkab(s1,s2,alpha,beta,i,j,k)[1]*spins::Sk(s2,beta,i,j,k)[0]);
                                }
                            }

                        }
                    }

                }
            }
        }
    }
    void hcpuConvFourier()
    {
        fields::hHk.IFill(0);

        //perform convolution in fourier space
        register unsigned int i = 0,j = 0, k = 0, alpha=0, beta=0;
        //convolute for each atomic plane
        for(i = 0 ; i < geom::dim[0]*geom::Nk[0] ; i++)
        {
            for(j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    for(alpha = 0 ; alpha < 3 ; alpha++)
                    {
                        for(beta = 0 ; beta < 3 ; beta++)
                        {
                            fields::hHk(alpha,i,j,k)[0]+=(intmat::hNkab(alpha,beta,i,j,k)[0]*spins::hSk(beta,i,j,k)[0]-intmat::hNkab(alpha,beta,i,j,k)[1]*spins::hSk(beta,i,j,k)[1]);
                            fields::hHk(alpha,i,j,k)[1]+=(intmat::hNkab(alpha,beta,i,j,k)[0]*spins::hSk(beta,i,j,k)[1]+intmat::hNkab(alpha,beta,i,j,k)[1]*spins::hSk(beta,i,j,k)[0]);
                            /*if(alpha==beta)
                            {
                                std::cout << "spins\t" << spins::hSk(beta,i,j,k)[0] << "\t" << spins::hSk(beta,i,j,k)[1] << std::endl;
                                std::cout << "IM\t" << intmat::hNkab(alpha,beta,i,j,k)[0] << "\t" << intmat::hNkab(alpha,beta,i,j,k)[1] << std::endl;

                    std::cout << fields::hHk(alpha,i,j,k)[0] << std::endl;
                            }*/
                        }
                    }

                }
            }
        }
    }
    void dipcpuConvFourier()
    {
        fields::dipHk.IFill(0);

        //perform convolution in fourier space
        register unsigned int i = 0,j = 0, k = 0, s1 = 0, s2=0, alpha=0, beta=0;
        for(i = 0 ; i < geom::zpdim[0]*geom::Nk[0] ; i++)
        {
            for(j = 0 ; j < geom::zpdim[1]*geom::Nk[1] ; j++)
            {
                for(k = 0 ; k < geom::zpdim[2]*geom::Nk[2] ; k++)
                {
                    for(alpha = 0 ; alpha < 3 ; alpha++)
                    {
                        for(beta = 0 ; beta < 3 ; beta++)
                        {
                            fields::dipHk(alpha,i,j,k)[0]+=(intmat::dipNkab(alpha,beta,i,j,k)[0]*spins::dipSk(beta,i,j,k)[0]-intmat::dipNkab(alpha,beta,i,j,k)[1]*spins::dipSk(beta,i,j,k)[1]);
                            fields::dipHk(alpha,i,j,k)[1]+=(intmat::dipNkab(alpha,beta,i,j,k)[0]*spins::dipSk(beta,i,j,k)[1]+intmat::dipNkab(alpha,beta,i,j,k)[1]*spins::dipSk(beta,i,j,k)[0]);
        //                    std::cout << fields::dipHk(alpha,i,j,k)[0] << "\t" << intmat::dipNkab(alpha,beta,i,j,k)[0] << "\t" << intmat::dipNkab(alpha,beta,i,j,k)[1] << "\t" << spins::dipSk(beta,i,j,k)[0] << "\t" << spins::dipSk(beta,i,j,k)[1] << std::endl;
        //                    std::cin.get();
                        }
                    }


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
    void fillfloat(int size,double* da,float* fa)
    {
        for(int i = 0 ; i < size ; i++)
        {
            fa[i]=float(da[i]);
        }
    }
    void fillfloat(int size0,int size1,int size2,Array3D<fftw_complex> da,Array3D<fftwf_complex> fa)
    {
        for(int i = 0 ; i < size0 ; i++)
        {
            for(int j = 0 ; j < size1 ; j++)
            {
                for(int k = 0 ; k < size2 ; k++)
                {
                    fa(i,j,k)[0]=static_cast<float>(da(i,j,k)[0]);
                    fa(i,j,k)[1]=static_cast<float>(da(i,j,k)[1]);
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
    void outputDiscVTU(unsigned int t)
    {
        std::stringstream pvss;
        std::string pv="discmap";
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
        pvf << "<Piece NumberOfPoints=\""<<(maxcx+1)*(maxcy+1)*(maxcz+1)<<"\"  NumberOfCells=\"1\">" << "\n";
        pvf << "<PointData Scalar=\"Spin\">" << "\n";
        pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        for(unsigned int i = 0 ; i < maxcx+1 ; i++)
        {
            for(unsigned int j = 0 ; j < maxcy+1 ; j++)
            {
                for(unsigned int k = 0 ; k < maxcz+1 ; k++)
                {
                    pvf << magdisc(i,j,k,0) << "\t" << magdisc(i,j,k,1) << "\t" << magdisc(i,j,k,2) << "\n";
                }
            }
        }
        pvf << "</DataArray>" << "\n";
        pvf << "</PointData>" << "\n";
        pvf << "<CellData>" << "\n";
        pvf << "</CellData>" << "\n";
        pvf << "<Points>" << "\n";
        pvf << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        for(unsigned int i = 0 ; i < maxcx+1 ; i++)
        {
            for(unsigned int j = 0 ; j < maxcy+1 ; j++)
            {
                for(unsigned int k = 0 ; k < maxcz+1 ; k++)
                {
                    pvf << static_cast<double>(i) << "\t" << static_cast<double>(j) << "\t" << static_cast<double>(k) << std::endl;
                }
            }
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

        if(llg::trans==0)
        {
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
        }
        else
        {
            pvf << "<?xml version=\"1.0\"?>" << "\n";
            pvf << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
            pvf << "<UnstructuredGrid>" << "\n";
            pvf << "<Piece NumberOfPoints=\""<<geom::dim[0]*geom::dim[1]*geom::dim[2]<<"\"  NumberOfCells=\"1\">" << "\n";
            pvf << "<PointData Scalar=\"Spin\">" << "\n";
            pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
            for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                    {
                        double op[3]={0,0,0};
                        for(unsigned int t = 0 ; t < geom::ucm.GetNMS() ; t++)
                        {
                            unsigned int atid=geom::atnolu(i,j,k,t);
                            op[0]+=spins::Sx[atid]*llg::optrans(t,0);
                            op[1]+=spins::Sy[atid]*llg::optrans(t,1);
                            op[2]+=spins::Sz[atid]*llg::optrans(t,2);
                        }
                        pvf << op[0] << "\t" << op[1] << "\t" << op[2] << std::endl;
                    }
                }
            }
        }

        pvf << "</DataArray>" << "\n";
        pvf << "</PointData>" << "\n";
        pvf << "<CellData>" << "\n";
        pvf << "</CellData>" << "\n";
        pvf << "<Points>" << "\n";
        pvf << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        if(llg::trans==0)
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                pvf << geom::rx[i] << "\t" << geom::ry[i] << "\t" << geom::rz[i] << "\n";
            }
        }
        else
        {
            for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                    {
                        pvf << static_cast<double>(i)*geom::rprim(0,0)+static_cast<double>(j)*geom::rprim(1,0)+static_cast<double>(k)*geom::rprim(2,0) << "\t";
                        pvf << static_cast<double>(i)*geom::rprim(0,1)+static_cast<double>(j)*geom::rprim(1,1)+static_cast<double>(k)*geom::rprim(2,1) << "\t";
                        pvf << static_cast<double>(i)*geom::rprim(0,2)+static_cast<double>(j)*geom::rprim(1,2)+static_cast<double>(k)*geom::rprim(2,2) << std::endl;
                    }
                }
            }
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
    void outputSpinsVTU(unsigned int t,double T)
    {
        std::stringstream pvss;
        std::string pv="spinmap";
        pvss << pv << "_T_" << T << "_t_" << t << ".vtu";
        std::string pvs = pvss.str();
        std::ofstream pvf(pvs.c_str());
        if(!pvf.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing paraview files");
        }
        if(llg::trans==0)
        {
            pvf << "<?xml version=\"1.0\"?>" << "\n";
            pvf << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
            pvf << "<UnstructuredGrid>" << "\n";
            pvf << "<Piece NumberOfPoints=\""<<geom::nspins<<"\"  NumberOfCells=\"1\">" << "\n";
            pvf << "<PointData Scalar=\"Spin\">" << "\n";
            pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";

            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                pvf << spins::Sx(i)*geom::mu[i] << "\t" << spins::Sy(i)*geom::mu[i] << "\t" << spins::Sz(i)*geom::mu[i] << "\n";
            }
        }
        else
        {
            pvf << "<?xml version=\"1.0\"?>" << "\n";
            pvf << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
            pvf << "<UnstructuredGrid>" << "\n";
            pvf << "<Piece NumberOfPoints=\""<<geom::dim[0]*geom::dim[1]*geom::dim[2]<<"\"  NumberOfCells=\"1\">" << "\n";
            pvf << "<PointData Scalar=\"Spin\">" << "\n";
            pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
            for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                    {
                        double op[3]={0,0,0};
                        for(unsigned int t = 0 ; t < geom::ucm.GetNMS() ; t++)
                        {
                            unsigned int atid=geom::atnolu(i,j,k,t);
                            op[0]+=spins::Sx[atid]*llg::optrans(t,0);
                            op[1]+=spins::Sy[atid]*llg::optrans(t,1);
                            op[2]+=spins::Sz[atid]*llg::optrans(t,2);
                        }
                        pvf << op[0] << "\t" << op[1] << "\t" << op[2] << std::endl;
                    }
                }
            }
        }

        pvf << "</DataArray>" << "\n";
        pvf << "</PointData>" << "\n";
        pvf << "<CellData>" << "\n";
        pvf << "</CellData>" << "\n";
        pvf << "<Points>" << "\n";
        pvf << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        if(llg::trans==0)
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                pvf << geom::rx[i] << "\t" << geom::ry[i] << "\t" << geom::rz[i] << "\n";
            }
        }
        else
        {
            for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                    {
                        pvf << static_cast<double>(i)*geom::rprim(0,0)+static_cast<double>(j)*geom::rprim(1,0)+static_cast<double>(k)*geom::rprim(2,0) << "\t";
                        pvf << static_cast<double>(i)*geom::rprim(0,1)+static_cast<double>(j)*geom::rprim(1,1)+static_cast<double>(k)*geom::rprim(2,1) << "\t";
                        pvf << static_cast<double>(i)*geom::rprim(0,2)+static_cast<double>(j)*geom::rprim(1,2)+static_cast<double>(k)*geom::rprim(2,2) << std::endl;
                    }
                }
            }
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
    void outputSpinsVTU(double T)
    {
        std::stringstream pvss;
        std::string pv="spinmap";
        pvss << pv << "_T_" << T << ".vtu";
        std::string pvs = pvss.str();
        std::ofstream pvf(pvs.c_str());
        if(!pvf.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing paraview files");
        }
        if(llg::trans==0)
        {
            pvf << "<?xml version=\"1.0\"?>" << "\n";
            pvf << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
            pvf << "<UnstructuredGrid>" << "\n";
            pvf << "<Piece NumberOfPoints=\""<<geom::nspins<<"\"  NumberOfCells=\"1\">" << "\n";
            pvf << "<PointData Scalar=\"Spin\">" << "\n";
            pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";

            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
//                pvf << spins::Sx(i) << "\t" << spins::Sy(i) << "\t" << spins::Sz(i) << "\n";
                pvf << spins::Sx(i)*geom::mu[i] << "\t" << spins::Sy(i)*geom::mu[i] << "\t" << spins::Sz(i)*geom::mu[i] << "\n";
            }
        }
        else
        {
            pvf << "<?xml version=\"1.0\"?>" << "\n";
            pvf << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
            pvf << "<UnstructuredGrid>" << "\n";
            pvf << "<Piece NumberOfPoints=\""<<geom::dim[0]*geom::dim[1]*geom::dim[2]<<"\"  NumberOfCells=\"1\">" << "\n";
            pvf << "<PointData Scalar=\"Spin\">" << "\n";
            pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
            for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                    {
                        double op[3]={0,0,0};
                        for(unsigned int t = 0 ; t < geom::ucm.GetNMS() ; t++)
                        {
                            unsigned int atid=geom::atnolu(i,j,k,t);
                            op[0]+=spins::Sx[atid]*llg::optrans(t,0);
                            op[1]+=spins::Sy[atid]*llg::optrans(t,1);
                            op[2]+=spins::Sz[atid]*llg::optrans(t,2);
                        }
                        pvf << op[0] << "\t" << op[1] << "\t" << op[2] << std::endl;
                    }
                }
            }
        }

        pvf << "</DataArray>" << "\n";
        pvf << "</PointData>" << "\n";
        pvf << "<CellData>" << "\n";
        pvf << "</CellData>" << "\n";
        pvf << "<Points>" << "\n";
        pvf << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        if(llg::trans==0)
        {
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {
                pvf << geom::rx[i] << "\t" << geom::ry[i] << "\t" << geom::rz[i] << "\n";
            }
        }
        else
        {
            for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
            {
                for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
                {
                    for(unsigned int k = 0 ; k < geom::dim[2] ; k++)
                    {
                        pvf << static_cast<double>(i)*geom::rprim(0,0)+static_cast<double>(j)*geom::rprim(1,0)+static_cast<double>(k)*geom::rprim(2,0) << "\t";
                        pvf << static_cast<double>(i)*geom::rprim(0,1)+static_cast<double>(j)*geom::rprim(1,1)+static_cast<double>(k)*geom::rprim(2,1) << "\t";
                        pvf << static_cast<double>(i)*geom::rprim(0,2)+static_cast<double>(j)*geom::rprim(1,2)+static_cast<double>(k)*geom::rprim(2,2) << std::endl;
                    }
                }
            }
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
    void calc_Ts()
    {
        for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
        {
            llg::Ts(s)=0.0;
            llg::cps(s)=0.0;
            llg::dps(s)=0.0;
        }
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            const double s[3]={spins::Sx[i],spins::Sy[i],spins::Sz[i]};
            const double h[3]={fields::Hx[i],fields::Hy[i],fields::Hz[i]};
            const double sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const double sdh=s[0]*h[0]+s[1]*h[1]+s[2]*h[2];
            //get the species
            unsigned int spec=geom::lu(i,3);
            double mm=geom::ucm.GetMu(spec);
            llg::cps(spec)+=mm*sxh[0]*sxh[0]+mm*sxh[1]*sxh[1]+mm*sxh[2]*sxh[2];
            llg::dps(spec)+=sdh;
        }
        for(unsigned int s = 0 ; s < geom::ucm.GetNMS() ; s++)
        {

            llg::Ts(s)=llg::muB*llg::cps(s)/(2.0*llg::kB*llg::dps(s));
            //std::cout << geom::ucm.GetMu(s) << std::endl;
        }
    }


}
