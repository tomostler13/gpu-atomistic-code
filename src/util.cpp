// File: util.cpp
// Author:Tom Ostler
// Last-modified: 04 Jan 2013 14:49:15
#include "../inc/util.h"
#include "../inc/spins.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/error.h"
#include "../inc/array.h"
#include "../inc/fields.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <sstream>
#include <fftw3.h>
//-------------------------------------
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace util
{
    bool ui=false;
    bool pvout=false;
    std::string pvf;
    unsigned int update=0;
    std::string dir;

    void initUtil(int argc,char *argv[])
    {
        assert(geom::gi);
        assert(fields::fi);
		config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Util details***" << std::endl;
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

		libconfig::Setting &setting = config::cfg.lookup("util");
        try
        {
		    setting.lookupValue("PVout",pvout);
        }
        catch(const libconfig::SettingTypeException &stex)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting type error");
        }

		FIXOUT(config::Info,"Output to paraview:" << config::isTF(pvout) << std::endl);
        try
        {
            if(pvout)
            {
                setting.lookupValue("file",pvf);
                setting.lookupValue("dir",dir);
                std::stringstream sstrdir;
                sstrdir << "mkdir " << dir;
                std::string str=sstrdir.str();
                system(str.c_str());
                setting.lookupValue("update",update);
            }
        }
        catch(const libconfig::SettingTypeException &stex)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting type error");
        }
        FIXOUT(config::Info,"Outputting to files:" << pvf << std::endl);
        FIXOUT(config::Info,"Update period:" << update*fields::dfu << " timesteps" << std::endl);
        ui=true;
    }

    void outputSpinsVTU(unsigned int t)
    {
        assert(geom::gi);
        assert(spins::si);
        assert(ui);
        std::stringstream pvss;
        pvss << dir << "/" << pvf << "_" << t << ".vtu";
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
        pvf << "<Piece NumberOfPoints=\""<<geom::ss<<"\"  NumberOfCells=\"1\">" << "\n";
        pvf << "<PointData Scalar=\"Spin\">" << "\n";
        pvf << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
                    pvf << spins::sx(i) << "\t" << spins::sy(i) << "\t" << spins::sz(i) << "\n";
        }
        pvf << "</DataArray>" << "\n";
        pvf << "</PointData>" << "\n";
        pvf << "<CellData>" << "\n";
        pvf << "</CellData>" << "\n";
        pvf << "<Points>" << "\n";
        pvf << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            pvf << double(geom::lu(i,0))*geom::gs[0] << "\t" << double(geom::lu(i,1))*geom::gs[1] << "\t" << double(geom::lu(i,2))*geom::gs[2] << "\n";
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
    std::string getdir(unsigned int i)
    {
        if(i==0)
        {
            return("x");
        }
        else if(i==1)
        {
            return("y");
        }
        else if(i==2)
        {
            return("z");
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Too many directions");
        }
    }

    void calcm()
    {
        spins::mx=0;
        spins::my=0;
        spins::mz=0;
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            spins::mx+=spins::sx[i];
            spins::my+=spins::sy[i];
            spins::mz+=spins::sz[i];
        }
        spins::modm=sqrt(spins::mx*spins::mx+spins::my*spins::my+spins::mz*spins::mz);
        spins::mx/=double(geom::ss);
        spins::my/=double(geom::ss);
        spins::mz/=double(geom::ss);
        spins::modm/=double(geom::ss);
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
}
