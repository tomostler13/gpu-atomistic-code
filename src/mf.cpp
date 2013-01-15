#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>
#include "../inc/array.h"
#include "../inc/array2d.h"
#include "../inc/mf.h"
#include "../inc/mat.h"
#include "../inc/config.h"
#include "../inc/tdp.h"
#include "../inc/geom.h"
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace mf
{
    bool mfi=false;
    Array2D<double> J0;
    //the identity matrix
    Array2D<double> I;
    double z=6;
    Array<double> fmsub;
    Array<double> D;
    Array<double> beta;
    Array<double> hmf;
    Array2D<double> hmfg;


    void initmf(int argc,char *argv[])
    {
        assert(mat::mi);
        config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Mean Field details***" << std::endl;
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

        libconfig::Setting &setting = config::cfg.lookup("mf");

        setting.lookupValue("mfz",z);
        J0.resize(mat::nlat,mat::nlat);
        I.resize(mat::nlat,mat::nlat);
        config::Info << "Mean field exchange interaction matrix details" << std::endl;
        for(unsigned int i = 0 ; i < mat::nlat ; i++)
        {
            std::stringstream sstr;
            sstr << "Tc" << i << "j";
            std::string str=sstr.str();
            config::Info << "J " << i << "\t\t\t[";

            for(unsigned int j = 0 ; j < mat::nlat ; j++)
            {
                J0(i,j)=setting[str.c_str()][j];
                //use the usual MF expression k_B T_c = J/3
                //and find the exchange constant
                J0(i,j)*=(3.0*1.38e-16);
                if(j < (mat::nlat-1))
                {
                    config::Info << J0(i,j) << ",";
                }
                else if (j==(mat::nlat-1))
                {
                    config::Info << J0(i,j) << "]" << std::endl;
                }

                if(i==j)
                {
                    I(i,j)=1.0;
                }
                else
                {
                    I(i,j)=0.0;
                }

            }

        }

        hmf.resize(mat::nlat);
        beta.resize(geom::ss);
        hmfg.resize(geom::ss,mat::nlat);
        D.resize(mat::nlat);
        fmsub.resize(mat::nlat);
        D.IFill(0);
        config::Info << "Mean Field anisotropy constants\nSublattice:" << std::endl;
        for(unsigned int i = 0 ; i < mat::nlat ; i++)
        {
            D[i]=setting["D"][i];
            D[i]*=1e-16;
            FIXOUT(config::Info,i << D[i] << std::endl);
        }
        fmsub.IFill(0);
        mfi=true;
    }
    void usrfun(double *msub,int n,double *fvec,double **fjac,double& lbeta,unsigned int& atom)
    {

        //calculate the sublattice magnetization Langevin function
        for(unsigned int i = 0 ; i < n ; i++)
        {
            hmf[i]=D[i];
            for(unsigned int j = 0 ; j < n ; j++)
            {
                hmf[i] += mat::conc[j]*J0(i,j)*msub[j+1];
            }
            hmfg(atom,i)=hmf[i];
            fvec[i+1]=msub[i+1]-langevin(hmf[i]*lbeta);
        }
        //calculate the Jacobian matrix
        for(unsigned int i = 0 ; i < n ; i++)
        {
            for(unsigned int j = 0 ; j < n ; j++)
            {
                fjac[i+1][j+1] = I(i,j);
                fjac[i+1][j+1] -= (  (2.0*D[i]*I(i,j)+mat::conc[j]*J0(i,j))*lbeta*dlangevin(hmf[i]*lbeta));
//                std::cout << lbeta << "\t" << hmf[i]*lbeta << "\t" << dlangevin(hmf[i]*lbeta) << std::endl;
//                std::cout << fjac[i+1][j+1] << std::endl;
            }
//            std::cin.get();
        }

    }

    double langevin(double x)
    {
        return((1./tanh(x)-(1./x)));
    }
    //derivative of the langevin function
    double dlangevin(double x)
    {
        return((1./(x*x))-(1./(sinh(x)*sinh(x))));
    }
}
