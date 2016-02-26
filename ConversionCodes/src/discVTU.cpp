// File: discVTU.cpp
// Author:Tom Ostler
// Created: 25 Feb 2016
// Last-modified: 25 Feb 2016 23:19:55

//The purpose of this section of code is to take the output of the vtu
//files and put the spins into cells and calculate the average magnetization
//within the cell. This aids in visualising systems with too many spins.
//At the moment it only works for a single species because as of the creation
//date the VTU output format of the spin dynamics code does not have any
//knowledge of the species
#include <cmath>
#include <iostream>
#include <libconfig.h++>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "../../inc/error.h"
#include "../../inc/array.h"
#include "../../inc/array2d.h"
#include "../../inc/array3d.h"
#include "../../inc/array4d.h"
#include "../../inc/defines.h"
#include "../../inc/random.h"
#include "../inc/inputs.h"
int main(int argc,char *argv[])
{
    Array3D<unsigned int> count;
    Array4D<double> mag;
    Array2D<double> spins;
    std::cout << "#Reading config file..." << std::flush;
    inputs::readcff(argc,argv);
    std::cout << "Done" << std::endl;

    std::cout << "#Resizing magnetization and spin arrays..." << std::flush;
    int ddim[3]={inputs::dimin[0]*inputs::Nk[0]/inputs::dimout[0],inputs::dimin[1]*inputs::Nk[1]/inputs::dimout[1],inputs::dimin[2]*inputs::Nk[2]/inputs::dimout[2]};
    mag.resize(ddim[0],ddim[1],ddim[2],3);
    mag.IFill(0);
    count.resize(ddim[0],ddim[1],ddim[2]);
    spins.resize(inputs::nspins,6);
    spins.IFill(0);
    std::cout << "Done" << std::endl;
    for(unsigned int t = inputs::lt ; t < inputs::ut+1 ; t+=inputs::dt)
    {
        std::stringstream sstr;
        sstr << inputs::fb << t << ".vtu";
        std::string str=sstr.str();

        std::cout << "#Opening file " << str << "..." << std::flush;
        std::ifstream ifs(str.c_str());
        if(!ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open input file");
        }
        std::cout << "done" << std::endl;
        std::string dumpline;
        //get rid of the top six lines
        for(unsigned int i = 0 ; i < 6 ; i++)
        {
            std::getline(ifs,dumpline);
        }
        //get the spin vectors (sx,sy,sz)
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            ifs >> spins(i,0) >> spins(i,1) >> spins(i,2);
        }
        //std::cout << spins(inputs::nspins-1,0) << std::endl;
        //get rid of the middle six lines
        for(unsigned int i = 0 ; i < 7 ; i++)
        {
            std::getline(ifs,dumpline);
        //    std::cout << dumpline << std::endl;
        }
        //get the coordinates
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            ifs >> spins(i,3) >> spins(i,4) >> spins(i,5);
        //    std::cout << spins(i,3) << "\t" << spins(i,4) << "\t" << spins(i,5) << std::endl;
        }
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        ifs.close();
        if(ifs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not close the ifs file");
        }
        count.IFill(0);
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        for(unsigned int i = 0 ; i < inputs::nspins ; i++)
        {
            unsigned int cx = static_cast<unsigned int>(spins(i,3)/inputs::abc[0])*inputs::Nk[0];
            unsigned int cy = static_cast<unsigned int>(spins(i,4)/inputs::abc[1])*inputs::Nk[1];
            unsigned int cz = static_cast<unsigned int>(spins(i,5)/inputs::abc[2])*inputs::Nk[2];
            unsigned int dx=static_cast<unsigned int>((static_cast<double>(cx)/static_cast<double>(inputs::dimout[0]))+0.);
            unsigned int dy=static_cast<unsigned int>((static_cast<double>(cy)/static_cast<double>(inputs::dimout[1]))+0.);
            unsigned int dz=static_cast<unsigned int>((static_cast<double>(cz)/static_cast<double>(inputs::dimout[2]))+0.);
            //std::cout << cx << "\t" << cy << "\t" << cz << "\t" << dx << "\t" << dy << "\t" << dz << std::endl;
            count(dx,dy,dz)++;//=count(dx,dy,dz)+1;
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
            mag(dx,dy,dz,0)+=spins(i,0);
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
            mag(dx,dy,dz,1)+=spins(i,1);
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
            mag(dx,dy,dz,2)+=spins(i,2);
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        }
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        std::stringstream osstr;
        osstr << "disc_" << t << ".vtu";
        std::string ostr=osstr.str();
        std::ofstream ofs;
        ofs.open(ostr.c_str());
        if(!ofs.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Couldn't open an output file");
        }
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        ofs << "<?xml version=\"1.0\"?>" << "\n";
        ofs << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
        ofs << "<UnstructuredGrid>" << "\n";
        ofs << "<Piece NumberOfPoints=\""<<ddim[0]*ddim[1]*ddim[2]<<"\"  NumberOfCells=\"1\">" << "\n";
        ofs << "<PointData Scalar=\"Spin\">" << "\n";
        ofs << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
//        std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        for(unsigned int dx = 0 ; dx < ddim[0] ; dx++)
        {
            for(unsigned int dy = 0 ; dy < ddim[1] ; dy++)
            {
                for(unsigned int dz = 0 ; dz < ddim[2] ; dz++)
                {
                    for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                    {
                        mag(dx,dy,dz,xyz)/=(static_cast<double>(count(dx,dy,dz)));
                    }
                    ofs << mag(dx,dy,dz,0) << "\t" << mag(dx,dy,dz,1) << "\t" << mag(dx,dy,dz,2) << std::endl;
                }
            }
        }
        ofs << "</DataArray>" << "\n";
        ofs << "</PointData>" << "\n";
        ofs << "<CellData>" << "\n";
        ofs << "</CellData>" << "\n";
        ofs << "<Points>" << "\n";
        ofs << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
        for(unsigned int dx = 0 ; dx < ddim[0] ; dx++)
        {
            for(unsigned int dy = 0 ; dy < ddim[1] ; dy++)
            {
                for(unsigned int dz = 0 ; dz < ddim[2] ; dz++)
                {
                    ofs << dx << "\t" << dy << "\t" << dz << std::endl;
                }
            }
        }
        ofs << "</DataArray>" << "\n";
        ofs << "</Points>" << "\n";
        ofs << "<Cells>" << "\n";
        ofs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n";
        ofs << "1" << "\n";
        ofs << "</DataArray>" << "\n";
        ofs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n";
        ofs << "1" << "\n";
        ofs << "</DataArray>" << "\n";
        ofs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << "\n";
        ofs << "1" << "\n";
        ofs << "</DataArray>" << "\n";
        ofs << "</Cells>" << "\n";
        ofs << "</Piece>" << "\n";
        ofs << "</UnstructuredGrid>" << "\n";
        ofs << "</VTKFile>" << "\n";
        ofs.close();
    }

    return(EXIT_SUCCESS);
}
