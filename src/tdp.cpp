// File: tdp.cpp
// Author:Tom Ostler
// Last-modified: 31 Dec 2012 13:56:37
#include "../inc/tdp.h"
#include "../inc/array.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/mat.h"
#include "../inc/fields.h"
#include "../inc/nrutil.h"
#include "../inc/mf.h"
#include "../inc/nrf.h"
#include <ctype.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
//==========================================================
#include <iostream>
#include <libconfig.h++>
#include <cstdlib>
#include <cstring>
#include <cmath>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace tdp
{
    bool tdpi=false;
    std::string toh;
    //this variable must be set if the system
    //temperature is uniform
    double uniformtemp=0;
    Array<double> systemp;
    Array<double> syschipar;
    Array<double> syschiperp;
    Array<double> sysexchstiff;
    Array<double> sysme;
    Array<double> sysalphaperp;
    Array<double> sysalphapar;
    Array<double> sysW1pf;
    Array<double> sysW2pf;
    Array<double> sysSigma;
    std::string chiperptype;
    std::string chipartype;
    std::string exchstifftype;
    std::string metype;
    unsigned int chiperpfunc;
    unsigned int chiparfunc;
    unsigned int exchstifffunc;
    unsigned int mefunc;
    //tss = transverse susecptibilitiy scale.
    //The transverse susceptibility can be scaled to
    //give something like the correct anisotropy
    //though the form will remain that of FePt.
    bool tss=false;
    //scaling factor (perp suscep)
    double tssf=0.0;
    void (*Tcalc)();
    void (*chiperpfp)();
    void (*chiparfp)();
    void (*exchstifffp)();
    void (*mefp)();

    void inittdp(int argc, char *argv[])
    {
        assert(geom::gi);
        assert(fields::fi);
        config::printline(config::Info);
		config::Info.width(45);config::Info << std::right << "*" << "**Temperature dependence details***" << std::endl;
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

        libconfig::Setting &setting = config::cfg.lookup("tdp");
		try
		{
			setting.lookupValue("toh",toh);
		}
		catch(const libconfig::SettingTypeException &stex)
		{
			error::errPreamble(__FILE__,__LINE__);
			error::errMessage("Setting type error");
		}
        FIXOUT(config::Info,"Type of heating:" << toh << std::endl);
        systemp.resize(geom::ss);
        syschipar.resize(geom::ss);
        syschiperp.resize(geom::ss);
        sysexchstiff.resize(geom::ss);
        sysme.resize(geom::ss);
        sysalphaperp.resize(geom::ss);
        sysalphapar.resize(geom::ss);
        sysW1pf.resize(geom::ss);
        sysW2pf.resize(geom::ss);
        sysSigma.resize(geom::ss);
        sysSigma.IFill(0);
        sysW1pf.IFill(0);
        sysW2pf.IFill(0);
        systemp.IFill(0);
        syschipar.IFill(0);
        syschiperp.IFill(0);
        sysexchstiff.IFill(0);
        sysme.IFill(1);
        sysalphaperp.IFill(0);
        sysalphapar.IFill(0);


        Tcalc=Tcalc0u;
        try
        {
            setting.lookupValue("chiperp",chiperptype);
            setting.lookupValue("chipar",chipartype);
            setting.lookupValue("exchstiffnesstype",exchstifftype);
            setting.lookupValue("chiperpfunction",chiperpfunc);
            setting.lookupValue("chiperfunction",chiparfunc);
            setting.lookupValue("exchstiffnessfunction",exchstifffunc);
            setting.lookupValue("metype",metype);
            setting.lookupValue("mefunction",mefunc);
        }
        catch(const libconfig::SettingTypeException &stex)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Setting type erro");
        }
        FIXOUT(config::Info,"Perpendicular susecptibility calculation kind:" << chiperptype << std::endl);
        if(chiperptype=="function")
        {
            FIXOUT(config::Info,"Function selected:" << chiperpfunc << std::endl);
            setting.lookupValue("scaleTransSuscep",tss);
            FIXOUT(config::Info,"Scaling the transverse susceptibility?" << config::isTF(tss) << std::endl);
            if(tss==true)
            {
                tssf=mat::Ms/(2.0*mat::K*chiperpsingle(0.0));
                FIXOUT(config::Info,"Scaling factor:" << tssf << std::endl);
            }
            else
            {
                tssf=1.0;
            }

        }
        else if(chipartype=="low temperature approximation")
        {
            //Do nothing at the moment
        }
        else if(chipartype=="meanfield")
        {
            //Do nothing at the moment
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Perpendicular susceptibility calculation type not recognised");
        }
        FIXOUT(config::Info,"Parallel susecptibility calculation kind:" << chipartype << std::endl);
        if(chipartype=="function")
        {
            FIXOUT(config::Info,"Function selected:" << chiparfunc << std::endl);
        }
        else if(chipartype=="meanfield")
        {
            //Do nothing at the moment
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Parallel susceptibility calculation type not recognised");
        }
        FIXOUT(config::Info,"Exchange stiffness calculation kind:" << exchstifftype << std::endl);
        if(exchstifftype=="function")
        {
            FIXOUT(config::Info,"Function selected:" << exchstifffunc << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Exchange stiffness calculation type not recognised");
        }
        FIXOUT(config::Info,"Equilibrium magnetization calculation kind:" << metype << std::endl);
        if(metype=="function")
        {
            FIXOUT(config::Info,"Function selected:" << mefunc << std::endl);
        }
        else if(metype=="meanfield")
        {
            //Do nothing at the moment
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Equilibrium magnetization calculation type not recognised");
        }
        if((metype=="meanfield" && chipartype=="function") || (metype=="function" && chipartype=="meanfield") || (metype=="meanfield" && chiperptype=="function") || (metype=="function" && chiperptype=="meanfield") || (chipartype=="function" && chiperptype=="meanfield") || (chipartype=="meanfield" && chiperptype=="function") )
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("Mixing of input function type and mean field type detected. Are you sure you want to do this?");
        }


        std::string ts = mat::materialid;
        for(unsigned int i = 0 ; i < ts.size() ; i++)
        {
            ts[i]=std::tolower(ts[i]);
        }

        if(chiperptype=="function" && chiperpfunc==0 && ts!="fept")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("***\n\n ***You have chosen the FePt functions from PRB, Vol. 77, 184428 (2008) as your input function  ***\n ***for the perpendicular susceptibility (chiperp) but your material id is not specified as FePt");
        }
        if(chipartype=="function" && chiparfunc==0 && ts!="fept")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("***\n\n ***You have chosen the FePt functions from PRB, Vol. 77, 184428 (2008) as your input function  ***\n ***for the parallel susceptibility (chipar) but your material id is not specified as FePt      ");
        }
        if(exchstifffunc==0 && ts!="fept")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errWarning("***\n\n ***You have chosen the FePt functions from PRB, Vol. 77, 184428 (2008) as your input function  ***\n ***for the exchange stiffness but your material id is not specified as FePt                    ");
        }
        if(metype=="function" && mefunc==0 && ts!="fept")
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("***\n\n ***You have chosen the FePt functions from PRB, Vol. 77, 184428 (2008) as your input function  ***\n ***for the equilibrium magnetization (me) but your material id is not specified as FePt        ");
        }

        if(chiperptype=="function" && chiperpfunc==0)
        {
            chiperpfp=chiperp0;
        }
        else if(chiperptype=="meanfield")
        {
            chiperpfp=chiperpmf;
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Chi perp method not recognized");
        }
        if(chipartype=="function" && chiparfunc==0)
        {
            chiparfp=chipar0;
        }
        else if(metype=="meanfield")
        {
            chiparfp=chiparmf;
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Chi parallel method not recognized");
        }
        if(exchstifffunc==0)
        {
            exchstifffp=exchstiff0;
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Exchange stiffness method not recognized");
        }
        if(fields::checkexchange==false)
        {
            exchstifffp=exchstiffNULL;
        }
        if(metype=="function" && mefunc==0)
        {
            mefp=me0;
        }
        else if(metype=="meanfield")
        {
            mefp=memf;
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Equilibrium magnetization method not recognized");
        }
        tdpi=true;
    }
    void Tcalc0u()
    {
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            systemp[i]=uniformtemp;
        }
    }

    void chiperp0()
    {
        //Fit Parameter
        double a0 = 0.00211549427182711;
        double a1 = 0.110224660661792;
        double a2 = -0.855153260915204;
        double a3 = 3.42088365387997;
        double a4 = -7.84821585896818;
        double a5 = 10.3247035469514;
        double a6 = -6.85608273303224;
        double a7 = 0.797453198330591;
        double a8 = 1.53787854178089;
        double a9 = -0.627128148404525;

        double chi_CGS = 0.0;
        //double chi_SI  = 0.0;
        double chi     = 0.0;
        double PI= 3.1415926535897932384626433832795028841971693993;;

        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            double T=systemp[i];

            if(T<(1.065*mat::Tc)) chi_CGS = a0+ a1*pow(pow(((1.068*mat::Tc-T)/(1.068*mat::Tc)),0.5),2.)+ a2*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),2.)+ a3*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),3.)+ a4*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),4.) + a5*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),5.) + a6*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),6.) + a7*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),7.)+ a8*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),8.) + a9*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),9.);
            else chi_CGS = (0.8*1.4/660.*mat::Tc)/(4*PI)/(T-mat::Tc);

            //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
            //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
            chi = chi_CGS*9.54393845712027; // (Tesla)

            syschiperp[i]=chi*tssf;
        }
    }

    void chiperpmf()
    {
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            if(systemp[i]<mat::Tc)
            {
                const double oome=1./tdp::sysme[i];
                syschiperp[i]=mat::chippf*oome*oome;
            }
            else
            {
                syschiperp[i]=syschipar[i];
            }
        }
    }
    double chiperpsingle(double T)
    {
        //Fit Parameter
        double a0 = 0.00211549427182711;
        double a1 = 0.110224660661792;
        double a2 = -0.855153260915204;
        double a3 = 3.42088365387997;
        double a4 = -7.84821585896818;
        double a5 = 10.3247035469514;
        double a6 = -6.85608273303224;
        double a7 = 0.797453198330591;
        double a8 = 1.53787854178089;
        double a9 = -0.627128148404525;

        double chi_CGS = 0.0;
        //double chi_SI  = 0.0;
        double chi     = 0.0;
        double PI= 3.1415926535897932384626433832795028841971693993;;

        if(T<(1.065*mat::Tc)) chi_CGS = a0+ a1*pow(pow(((1.068*mat::Tc-T)/(1.068*mat::Tc)),0.5),2.)+ a2*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),2.)+ a3*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),3.)+ a4*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),4.) + a5*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),5.) + a6*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),6.) + a7*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),7.)+ a8*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),8.) + a9*pow((((1.068*mat::Tc)-T)/(1.068*mat::Tc)),9.);
        else chi_CGS = (0.8*1.4/660.*mat::Tc)/(4*PI)/(T-mat::Tc);

        //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
        //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
        chi = chi_CGS*9.54393845712027; // (Tesla)
        return(chi);
    }
    void chipar0()
    {
        //Fit Parameter
        double a0 = 0.8;
        double a1 =-2.2e-07;
        double a2 = 1.95e-13;
        double a3 =-1.3e-17;
        double a4 =-4e-23;
        double a5 =-6.5076312364e-32;

        double chi_CGS = 0.0;
        //double chi_SI  = 0.0;
        double chi = 0.0;
        double PI= 3.1415926535897932384626433832795028841971693993;;
        //double factor =  0.75947907;

        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            double T=systemp[i];
            double factor=mat::Tc-T;

            if(T<mat::Tc) chi_CGS =(a0/660.*mat::Tc)/(4.*PI)/(factor)+a1*factor+ a2*factor*factor*factor+a3*factor*factor*factor*factor+ a4*factor*factor*factor*factor*factor*factor+ a5*factor*factor*factor*factor*factor*factor*factor*factor*factor;
            else chi_CGS = (1.1*1.4/660.*mat::Tc)/(4*PI)/(T-mat::Tc);
            //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
            //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
            chi = chi_CGS*9.54393845712027+0.308e-14; // (Tesla)
            syschipar[i]=chi;
        }
    }
    //when the parallel susceptibility is evaluated in
    //the mean field the magnetization must first be
    //calculated.
    void chiparmf()
    {
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            const double lb=mf::beta[i];
            const double Bp=mf::dlangevin(mf::hmfg(i,0)*lb);
            const double J=mf::J0(0,0);
            //There is a factor of 10^4 here because
            //1T=10^4 Oe
            syschipar[i]=(mat::mu[0]*9.27e-21*1e4/J)*Bp*lb*J/(1.0-Bp*lb*J);
        }

    }
    void exchstiff0()
    {
        const double ex0 = 3.90858143659231e-13;
        const double ex1 = 5.65571902911896e-11;
        const double ex2 = -1.11221431025254e-10;
        const double ex3 = 1.67761522644194e-10;
        const double ex4 = -1.38437771856782e-10;
        const double ex5 = 4.6483423884759e-11;

        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            const double T=systemp[i];
            const double factor=(mat::Tc-T)/mat::Tc;
            double A=0.0;

            if(T<mat::Tc)
            {
                A=ex0+ex1*factor+ex2*factor*factor+ex3*factor*factor*factor+ex4*factor*factor*factor*factor+ex5*factor*factor*factor*factor*factor;
            }
            sysexchstiff[i]=A;
        }

    }
    void exchstiffNULL()
    {
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            sysexchstiff[i]=0;
        }
    }
    void me0()
    {
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            const double T=systemp[i];
            const double TCmT=mat::Tc-T;
            const double Tc_m_T_o_Tc=TCmT/mat::Tc;
            const double factor=(1.068*mat::Tc-T)/(1.068*mat::Tc);

            if(T<mat::Tc)
            {
                sysme[i] = 1.3*sqrt(Tc_m_T_o_Tc)-0.12*Tc_m_T_o_Tc-0.51*Tc_m_T_o_Tc*Tc_m_T_o_Tc + 0.37*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc - 0.01*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc-0.03*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc;
            }
            else
            {
                sysme[i]=0;
            }
        }
    }
    //mean field
    void memf()
    {

        double tolx=1e-5,tolf=1e-5;
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            double *x;
            x=dvector(1,mat::nlat);
            for(unsigned int j = 1 ; j <= mat::nlat ; j++)
            {
                //set the starting value of the trial magnetization equal to
                //the magnetization at the previous timestep
                x[j]=sysme[i];
            }
            //newton-rhapson method (from Numerical Recipes)
            int ntrial=1000;
            nrf::mnewt(ntrial,x,mat::nlat,tolx,tolf,i);
//            if(i==0)
//            {
//                std::cout << systemp[0] << "\t" << x[1] << "\t" << x[2] << std::endl;
//            }
            //This is because at this stage the mean field will only
            //be used for ferromagnets and we are only considering
            //the ferromagnetic Landau-Lifshitz-Bloch equation at the
            //moment
            sysme[i]=x[1];
            free_dvector(x,1,mat::nlat);
        }
    }

    void calcalphas()
    {
        for(unsigned int i = 0 ; i < geom::ss ; i++)
        {
            const double T=systemp[i];
            sysalphapar[i]=mat::lambda*2.0*T/(3.0*mat::Tc);

            if(T<=mat::Tc)
            {
                sysalphaperp[i]=mat::lambda*(1.0-(T/(3.0*mat::Tc)));
            }
            else
            {
                sysalphaperp[i]=sysalphapar[i];
            }

            sysW1pf[i]=sqrt(T*sysalphapar[i]);
            sysW2pf[i]=sqrt(T*(sysalphaperp[i]-sysalphapar[i]));
        }

    }


}
