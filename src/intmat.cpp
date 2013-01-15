// File: intmat.cpp
// Author:Tom Ostler
// Last-modified: 27 Dec 2012 11:53:12
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <libconfig.h++>
#include "../inc/intmat.h"
#include "../inc/array3d.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/mat.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace intmat
{
	Array3D<fftw_complex> Nxx;
	Array3D<fftw_complex> Nxy;
	Array3D<fftw_complex> Nxz;
	Array3D<fftw_complex> Nyx;
	Array3D<fftw_complex> Nyy;
	Array3D<fftw_complex> Nyz;
	Array3D<fftw_complex> Nzx;
	Array3D<fftw_complex> Nzy;
	Array3D<fftw_complex> Nzz;
    bool imi=false;
    bool imft=false;

	void initIntmat(int argc,char *argv[])
	{
		assert(geom::gi);
		config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Interaction matrix details***" << std::endl;
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

        Nxx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nxy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nxz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nyx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nyy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nyz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nzx.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nzy.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nzz.resize(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2]);
        Nxx.IFill(0.0);
        Nxy.IFill(0.0);
        Nxz.IFill(0.0);
        Nyx.IFill(0.0);
        Nyy.IFill(0.0);
        Nyz.IFill(0.0);
        Nzx.IFill(0.0);
        Nzy.IFill(0.0);
        Nzz.IFill(0.0);
        imi=true;

	}

	void fillIntmat()
	{
		FIXOUT(config::Info,"Filling interaction matrix" << std::flush);
		for(unsigned int i = 0 ; i < geom::zpdim[0] ; i++)
		{
			for(unsigned int j = 0 ; j < geom::zpdim[1] ; j++)
			{
				for(unsigned int k = 0 ; k < geom::zpdim[2] ; k++)
				{
					if(!(i==0 && j==0 && k==0) && i!=geom::dim[0] && j!=geom::dim[1] && k!=geom::dim[2])
					{
						int lc[3]={i,j,k};
						int tc[3]={i,j,k};
						//vector from i to j
						double rij[3]={0,0,0};
						for(unsigned int l = 0 ; l < 3 ; l++)
						{
							if(lc[l]>geom::dim[l])
							{
								lc[l]=geom::dim[l]-lc[l];
								tc[l]=geom::zpdim[l]+lc[l];
							}
							//find rij in m
							rij[l]=double(lc[l])*geom::gs[l];
						}
						double mrij=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
						double oomrij3=1./(mrij*mrij*mrij);
						//unit vector from i to j
						double eij[3]={rij[0]/mrij,rij[1]/mrij,rij[2]/mrij};
                        double lnxx=1e-7*((3.0*eij[0]*eij[0])-1.0)*oomrij3;
                        Nxx(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnxx*geom::gsV*mat::Ms;//(i,j,k)

                        double lnyy=1e-7*((3.0*eij[1]*eij[1])-1.0)*oomrij3;
                        Nyy(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnyy*geom::gsV*mat::Ms;//(i,j,k)

                        double lnzz=1e-7*((3.0*eij[2]*eij[2])-1.0)*oomrij3;
                        Nzz(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnzz*geom::gsV*mat::Ms;//(i,j,k)


                        //Now w_{ij}^{\Gamma \Lambda}=\frac{3e_{ij}^\Gamma e_{ij}^{\Lambda}}{r_{ij}^3}
                        double lnxy=1e-7*3.0*eij[1]*eij[0]*oomrij3;
                        Nxy(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnxy*geom::gsV*mat::Ms; //(i,j,k)

                        double lnxz=1e-7*3.0*eij[2]*eij[0]*oomrij3;
                        Nxz(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnxz*geom::gsV*mat::Ms; //(i,j,k)

                        double lnyx=1e-7*3.0*eij[0]*eij[1]*oomrij3;
                        Nyx(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnyx*geom::gsV*mat::Ms; //(i,j,k)

                        double lnyz=1e-7*3.0*eij[2]*eij[1]*oomrij3;
                        Nyz(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnyz*geom::gsV*mat::Ms; //(i,j,k)

                        double lnzx=1e-7*3.0*eij[0]*eij[2]*oomrij3;
                        Nzx(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnzx*geom::gsV*mat::Ms; //(i,j,k)

                        double lnzy=1e-7*3.0*eij[1]*eij[2]*oomrij3;
                        Nzy(tc[0]             ,tc[1]             ,tc[2]             )[0]=lnzy*geom::gsV*mat::Ms; //(i,j,k)

					}
				}
			}
		}
		config::Info << "Done" << std::endl;
	}

	void fftIntmat()
	{
        assert(imi);
		FIXOUT(config::Info,"Setting up fftw of interaction matrix:" << std::flush);
        //plans for the transform. Only done once, so not persistent
        fftw_plan NxxP,NxyP,NxzP,NyxP,NyyP,NyzP,NzxP,NzyP,NzzP;

        //the demag tensor 3d fftw plans
        NxxP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nxx.ptr(),Nxx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NxyP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nxy.ptr(),Nxy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NxzP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nxz.ptr(),Nxz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NyxP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nyx.ptr(),Nyx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NyyP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nyy.ptr(),Nyy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NyzP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nyz.ptr(),Nyz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NzxP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nzx.ptr(),Nzx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NzyP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nzy.ptr(),Nzy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
        NzzP=fftw_plan_dft_3d(geom::zpdim[0],geom::zpdim[1],geom::zpdim[2],Nzz.ptr(),Nzz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);


		config::Info << "Done" << std::endl;
		FIXOUT(config::Info,"Performing forward transform of N:" << std::flush);
                //execute the demag tensor transforms and destroy the plans
        fftw_execute(NxxP);
        fftw_execute(NxyP);
        fftw_execute(NxzP);
        fftw_execute(NyxP);
        fftw_execute(NyyP);
        fftw_execute(NyzP);
        fftw_execute(NzxP);
        fftw_execute(NzyP);
        fftw_execute(NzzP);

        fftw_destroy_plan(NxxP);
        fftw_destroy_plan(NxyP);
        fftw_destroy_plan(NxzP);
        fftw_destroy_plan(NyxP);
        fftw_destroy_plan(NyyP);
        fftw_destroy_plan(NyzP);
        fftw_destroy_plan(NzxP);
        fftw_destroy_plan(NzyP);
        fftw_destroy_plan(NzzP);
        imft=true;

		config::Info << "Done" << std::endl;
	}
}
