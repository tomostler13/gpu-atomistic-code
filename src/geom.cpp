// File: geom.cpp
// Author:Tom Ostler
// Last-modified: 15 Jan 2013 17:38:42
#include "../inc/geom.h"
#include "../inc/neigh.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/array3d.h"
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/util.h"
#include "../inc/ellipse_grid.hpp"
#include <libconfig.h++>
#include <sstream>
#include <string>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace geom
{
    //dimensions of system that exist (i.e. non-zero padded) in unit cells
    unsigned int dim[3]={0,0,0};
    //dimensions when zero padding is taken into account
    int zpdim[3]={0,0,0};
    //number of zero pad units
    int zps=0;
    //number of macrospins that exist
    unsigned int ss=0;
    //is geometry enabled?
	bool gi=false;
    std::string structure;
    //number of atoms in the unit cell
    unsigned int nauc=0;
    //positions of my atoms in my unit cell
    Array2D<double> atompos;
    //scaling factor to take the unit cell to integer coordinates
    Array<unsigned int> scalecoords;
    //number of spatial cells in total including those
    //that don't necessarily exist.
    unsigned int maxss=0;

    //When performing the numerical integration we loop over the number of atoms
    //in order. For the fourier transform method of calculating the dipolar
    //field we will copy to the spin array as this is probably more efficient for
    //the GPU code.
    //spin lookup. Given a spin number this returns the coordinates.
    Array2D<unsigned int> lu;
    //given an i,j,k this returns the atom number (as in lu).
    Array3D<int> coords;

    void initGeom(int argc,char *argv[])
    {

        assert(config::lcf);
        config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Geometry details***" << std::endl;
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

        libconfig::Setting &setting = config::cfg.lookup("system");
        setting.lookupValue("structure",structure);
        //output stream for writing structure
        std::ofstream outstruc;
        outstruc.open("outstruc.dat");
        if(!outstruc.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not open file for writing structure");
        }
        scalecoords.resize(3);
        for(unsigned int i = 0 ; i < 3 ; i++)
        {
            try
            {
                dim[i]=setting["dim"][i];
                scalecoords[i]=setting["scalecoords"][i];

            }
            catch(const libconfig::SettingTypeException &stex)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Setting type error");
            }
        }
        FIXOUT(config::Info,"Scale coords: " << "(" << scalecoords[0] << "," << scalecoords[1] << "," << scalecoords[2] << ")" << std::endl);
        //get the number of atoms in the unit cell
        nauc=setting["nauc"];
        atompos.resize(nauc,3);
        //read the atom positions in the unit cell
        for(unsigned int i = 0 ; i < nauc ; i++)
        {
            std::stringstream sstr;
            sstr << "atom" << i;
            std::string str=sstr.str();
            for(unsigned int j = 0 ; j < 3 ; j++)
            {
                atompos(i,j)=setting[str.c_str()][j];
            }
        }
        FIXOUT(config::Info,"Number of atoms in the unit cell:" << nauc << std::endl);
        for(unsigned int i = 0 ; i < nauc ; i++)
        {
            std::stringstream sstr;
            sstr << "Position vector for atom " << i << "(;
            std::string str=sstr.str();
            FIXOUT(config::Info,str.c_str() << atompos(i,0) << "," << atompos(i,1) << "\t" << atompos(i,2) << ")" << std::endl);
        }
        if(structure=="cuboid")
        {
            ss=dim[0]*dim[1]*dim[2]*nauc;
            maxss=ss;
        }
        else if(structure=="cylinder")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Cylinder not possible at the moment with atomistic spin dynamics. Some code exists but it needs checking.");
            //check that dim[0] and dim[1] are odd for cylinder
            if(dim[0]%2==0 || dim[1]%2==0 || dim[0] < 5 || dim[1] < 5 || dim[1]>dim[0])
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Dimensions for cylinder not correct, must be odd and 5 or greater in x and y dimension.\nAlso dim[0] should be >= dim[1]");
            }
            maxss=dim[0]*dim[1]*dim[2];
            //ss is not know yet.
        }
        //we are going to re-use the coords array.
        coords.resize(dim[0]*scalecoords[0],dim[1]*scalecoords[1],dim[2]*scalecoords[2]);
        for(unsigned int i = 0 ; i < dim[0]*scalecoords[0] ; i++)
        {
            for(unsigned int j = 0 ; j < dim[1]*scalecoords[1] ; j++)
            {
                for(unsigned int k = 0 ; k < dim[2]*scalecoords[2] ; k++)
                {
                    coords(i,j,k)=0;
                }
            }
        }
        //returns the points
        double *xy;
        if(structure=="cylinder")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Cylinder not possible at the moment with atomistic spin dynamics. Some code exists but it needs checking.");

            //This section of code uses the ELLIPSE_GRID routine from http://people.sc.fsu.edu/~jburkardt/cpp_src/ellipse_grid/ellipse_grid.html
            //radius of ellipse
            double r[2]={double(dim[0])/double(dim[1]),1.0};
            //centre of ellipse
            double c[2]={r[0],r[1]};
            //needs an integer giving the number of points on the minor half axis
            int n=int(dim[1]/2);//int(dim[0]/2);
            //number of grid points
            int ng=0;
            ng=ellipse_grid_count(n,r,c);
            ss=ng;
            xy=ellipse_grid(n,r,c,ng);
            //find the biggsest x and y val
            double big[2]={0,0};
            for(unsigned int i = 0 ; i < ss ; i++)
            {
                for(unsigned int j = 0 ; j < 2 ; j++)
                {
                    if(xy[2*i+j]>big[j]){big[j]=xy[2*i+j];}
                }
            }
            //scale the x and y coords
            for(unsigned int i = 0 ; i < ss ; i++)
            {
                for(unsigned int j = 0 ; j < 2 ; j++)
                {
                    xy[2*i+j]/=big[j];
                    xy[2*i+j]*=double(dim[j]);

                }

                for(unsigned int k = 0 ; k < dim[2] ; k++)
                {
                    coords(int(xy[2*i]+0.5)-1,int(xy[2*i+1]+0.5)-1,k)=1;
                }

            }

        }

        try
        {
            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                neigh::pbc[i] = setting["pbc"][i];
                if(neigh::pbc[i]==true)
                {
                    if(structure=="cylinder")
                    {
                        if(neigh::pbc[i]==true)
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("With cylinder structure you cannot have periodic boundary conditions");
                        }
                    }
                    if(geom::dim[i]==1)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream sstr;
                        sstr<<"It looks as though you are trying to have periodic boundaries with only 1 atom in the " << util::getdir(i) << " direction" << std::endl;
                        std::string str=sstr.str();
                        error::errWarning(str.c_str());
                    }
                }
            }
		}
		catch(const libconfig::SettingTypeException &stex)
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Setting type error");
		}
        FIXOUT(config::Info,"System structure:" << structure << std::endl);
        FIXOUT(config::Info,"System dimensions (all):" << "(" << dim[0] << "," << dim[1] << "," << dim[2] << ")" << std::endl);
        zpdim[0]=2*dim[0]*scalecoords[0];
        zpdim[1]=2*dim[1]*scalecoords[1];
        zpdim[2]=2*dim[2]*scalecoords[2];
        zps=zpdim[0]*zpdim[1]*zpdim[2];*scalecoords[0]*scalecoords[1]*scalecoords[2];
        //should be multiplied by dim[2] as ss so far is in 2D
        lu.resize(ss*dim[2],4);
        lu.IFill(0);
        Array3D<int> tc;
        tc.resize(dim[0]*scalecoords[0],dim[1]*scalecoords[1],dim[2]*scalecoords[2]);
        unsigned int counter=0;
        for(unsigned int i = 0 ; i < geom::dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < geom::dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < geom::dim[2] ;  k++)
                {
                    //this can be used for debugging the neighbour list quickly.
                    //Plot this in gnuplot using
                    //sp "<filename>" u 2:3:4:1 with labels
                    /*  FOR DEBUGGING
                    std::cout << "\"" << counter << "\"\t" << i << "\t" << j << "\t" << k << std::endl;
                    */
                    if(structure=="cuboid")
                    {
                        for(unsigned int uc = 0 ; uc < nauc ; uc++)
                        {
                            lu(counter,0)=(int(atompos(uc,0)*double(scalecoords[0])))*i;
                            lu(counter,1)=(int(atompos(uc,1)*double(scalecoords[1])))*j;
                            lu(counter,2)=(int(atompos(uc,2)*double(scalecoords[2])))*k;
                            //atom exists
                            lu(counter,3)=1;
                            coords(i,j,k)=counter;
                            outstruc << i << "\t" << j << "\t" << k << std::endl;
                            counter++;
                        }
                    }
                    else if(structure=="cylinder")
                    {
                        //pre-determined by ellipse_grid
                        if(coords(i,j,k)==1)
                        {
                            lu(counter,0)=i;
                            lu(counter,1)=j;
                            lu(counter,2)=k;
                            lu(counter,3)=1;
                            tc(i,j,k)=counter;
                            counter++;
                            outstruc << "\"" << counter << "\"\t" << i << "\t" << j << "\t" << k << std::endl;

                        }
                        else
                        {
                            tc(i,j,k)=-1;
                        }

                    }
                }
            }
        }
        if(structure=="cuboid" && maxss!=ss)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Something is wrong, maxss should equal ss for cuboid");
        }

        for(unsigned int i = 0 ; i < dim[0] ; i++)
        {
            for(unsigned int j = 0 ; j < dim[1] ; j++)
            {
                for(unsigned int k = 0 ; k < dim[2] ; k++)
                {
                    coords(i,j,k)=tc(i,j,k);
                }
            }
        }
        tc.clear();

        FIXOUT(config::Info,"Number of macrospins: " << ss << std::endl);
        FIXOUT(config::Info,"Zero pad size (for FFT):" << zps << std::endl);
        FIXOUT(config::Info,"Grain size dimensions:" << "(" << gs[0] << "," << gs[1] << "," << gs[2] << ") m" << std::endl);
        //set the grain size volumne
        gsV=gs[0]*gs[1]*gs[2];
        FIXOUT(config::Info,"Grain size volume:" << gsV << " m^3" << std::endl);

		assert(zpdim[0]>0);
		assert(zpdim[1]>0);
		assert(zpdim[2]>0);
        assert(dim[0]>0);
        assert(dim[1]>0);
        assert(dim[2]>0);
        assert(cplxzpdim>0);
		gi=true;
        outstruc.close();
        if(outstruc.is_open())
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Could not close structure file");
        }
    }
}
