// File: llg.cpp
// Author:Tom Ostler
// Created: 22 Jan 2013
// Last-modified: 19 May 2023 13:20:54
#include "../inc/llg.h"
#include "../inc/llgCPU.h"
#include "../inc/config.h"
#include "../inc/defines.h"
#include "../inc/geom.h"
#include "../inc/spins.h"
#include "../inc/util.h"
#include "../inc/rscf.h"
#include "../inc/fields.h"
#include <cmath>
#include <sstream>
#include <string>
#include <cctype>
#include <algorithm>
#ifdef CUDA
#include <cuda.h>
#include "../inc/cuda.h"
#endif /*CUDA*/
namespace llg
{
    double applied[3]={0,0,0},T,dt,rdt,gyro=1.76e11,muB=9.27e-24,kB=1.38e-23;
    Array2D<double> optrans;
    unsigned int trans=0,intscheme;
    bool exstr=false;
    Array<double> Ts,dps,cps;
    Array<double> llgpf;
    std::string scm,intschemestr;

    void initLLG()
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**LLG details***" << std::endl;

        bool llgset=false;
        if(!config::cfg.exists("llg"))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Setting class llg does not exists, setting defaults.");
            llgset=false;
            dt=0.1e-15;
            FIXOUT(config::Info,"Timestep (default):" << dt << " seconds" << std::endl);
            rdt=dt*gyro;
            FIXOUT(config::Info,"Reduced timestep (default):" << rdt << std::endl);
            applied[0]=0.0;applied[1]=0.0;applied[2]=0.0;
            FIXOUTVEC(config::Info,"Applied field (default):",applied[0],applied[1],applied[2]);
            spins::update=100;
            FIXOUT(config::Info,"Update of spins (default):" << spins::update << " [Timesteps]" << std::endl;);
            spins::mag_calc_method=0;
            FIXOUT(config::Info,"Magnetization calculatioin method (default):" << spins::mag_calc_method << std::endl);
            scm="file";
            spins::setSpinsConfig();
            FIXOUT(config::Info,"Initial spin configuration (default):" << scm << std::endl);
            FIXOUT(config::Info,"Order parameter transform (default):" << trans << std::endl);
            util::escheck=false;
            FIXOUT(config::Info,"Outputting exchange striction (default):" << config::isTF(util::escheck) << std::endl);
            FIXOUT(config::Info,"Numerical integration method (default):" << "Heun" << std::endl);

        }
        else
        {
            libconfig::Setting &setting = config::cfg.lookup("llg");
//            std::cout << "does it exist or not " << config::cfg.exists("llg.exch_striction") << std::endl;
//            std::cin.get();
            if(!config::cfg.exists("llg.exch_striction"))
            {
                util::escheck=false;
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Setting, llg.exch_striction not set, defaulting to false.");
                FIXOUT(config::Info,"Outputting exchange striction (default):" << config::isTF(util::escheck) << std::endl);
            }
            else
            {
                util::escheck=true;
                FIXOUT(config::Info,"Outputting exchange striction:" << config::isTF(util::escheck) << std::endl);
                //open the file to write the polarisation to
                util::outESP.open("exch_stric_P.dat");
                if(!util::outESP.is_open())
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Could not open file exch_stric_P.dat, check disk space and access");
                }
                else
                {
                    util::outESP << "# This file contains the time in seconds in the first column followed by the" << std::endl;
                    util::outESP << "# dot product between spin pairs (each pair is output) normalised by the number" << std::endl;
                    util::outESP << "# of unit cell (N) and the number of interactions that each spin participates in" << std::endl;
                }

                if(util::escheck)//then read the variables that control it
                {
                    libconfig::Setting &esset = config::cfg.lookup("llg.exch_striction");
                    //get the number of pairs
                    if(!esset.lookupValue("NumPairs",util::esnp))
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("You have not correctly set up the exchange striction calculation. llg.exch_striction.NumPairs (int) must be set.");
                    }
                    else
                    {
                        FIXOUT(config::Info,"Number of pairs in exchange striction term:" << util::esnp << std::endl);
                        if(!esset.lookupValue("MaxUnitCells",util::maxuc))
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("You must set the maximum number of unit cells to look for for the exchange striction. Setting llg.exch_striction.MaxUnitCells (int)");
                        }
                        else
                        {
                            FIXOUT(config::Info,"Max number of unit cells to look for for the exchange striction:" << util::maxuc << std::endl);
                        }
                        util::espairs.resize(util::esnp,2);
                        util::espairs.IFill(0);
                        //resize the arrays to hold the unit cell interactio details
                        util::esnum.resize(util::esnp);
                        util::esnum.IFill(0);
                        util::esuc.resize(util::esnp,util::maxuc,3);
                        util::esuc.IFill(0);
                        util::lambda.resize(util::esnp);
                        util::lambda.IFill(0);
                        util::esP.resize(util::esnp);
                        util::esP.IFill(0);

                    }
                    for(unsigned int i = 0 ; i < util::esnp ; i++)
                    {
                        std::stringstream sstr;
                        sstr << "Pair" << i;
                        std::string str=sstr.str();

                        for(unsigned int  j = 0 ; j < 2 ; j++)
                        {
                            try
                            {
                                util::espairs(i,j)=esset[str.c_str()][j];
                            }
                            catch(const libconfig::SettingNotFoundException &snf)
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                std::stringstream errsstr;
                                errsstr << "Setting not found exception caught. Exchange striction not set up properly. Setting " << snf.getPath() << " must be set.";
                                std::string errstr=errsstr.str();
                                error::errMessage(errstr);
                            }
                        }
                        std::stringstream opsstr,opsstr1;
                        opsstr << "Pair " << i << " is between atoms:";
                        opsstr1 << util::espairs(i,0) << " and " << util::espairs(i,1);
                        std::string opstr=opsstr.str(),opstr1=opsstr1.str();
                        FIXOUT(config::Info,opstr << opstr1 << std::endl);
                        opsstr.str("");
                        opsstr.clear();
                        opsstr << "llg.exch_striction.UnitCells" << i;
                        opstr=opsstr.str();
                        if(!config::cfg.exists(opstr.c_str()))
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            opsstr.str("");
                            opsstr.clear();
                            opsstr << "Your setting of the unit cell directions llg.exch_striction.UnitCell" << i << " cannout be found, check your config file.";
                            opstr=opsstr.str();
                            error::errMessage(opstr);
                        }
                        else
                        {
                            libconfig::Setting &pset = config::cfg.lookup(opstr.c_str());

                            if(!pset.lookupValue("lambda",util::lambda(i)))
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                opsstr.str("");opsstr.clear();
                                opsstr << "The exchange striction constant for pair " << i << " was not set properly, setting llg.exch_striction.UnitCells" << i << ".lambda (double)";
                                error::errMessage(opsstr);
                            }
                            else
                            {
                                FIXOUT(config::Info,"Exchange striction constant:" << util::lambda(i) << " [UNITS????]" << std::endl);
                            }
                            if(!pset.lookupValue("num",util::esnum(i)))
                            {
                                opsstr.str("");
                                opsstr.clear();
                                opsstr << "Could not read number of unit cells for calculation of exchange striction (setting llg.exch_striction.UnitCell" << i << ".num (int)";
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage(opsstr);
                            }
                            else
                            {

                                FIXOUT(config::Info,"Number of unit cells:" << util::esnum(i) << std::endl);
                                if(util::esnum(i)>util::maxuc)
                                {
                                    error::errPreamble(__FILE__,__LINE__);
                                    opsstr.str("");opsstr.clear();
                                    opsstr << "The number of unit cells to look for (llg.exch_striction.UnitCells." << i << ".num (int) cannot be larger than setting llg.exch_striction.NumPairs (int), check your config file.";
                                    error::errMessage(opsstr);
                                }
                                for(unsigned int uc = 0 ; uc < util::esnum(i) ; uc++)
                                {
                                    opsstr.str("");opsstr.clear();
                                    opsstr << "dir" << uc;
                                    opstr=opsstr.str();
                                    for(unsigned int d = 0 ; d < 3 ; d++)
                                    {
                                        try
                                        {
                                            util::esuc(i,uc,d)=pset[opstr.c_str()][d];
                                        }
                                        catch(const libconfig::SettingNotFoundException &snf)
                                        {
                                            error::errPreamble(__FILE__,__LINE__);
                                            std::stringstream errsstr;
                                            errsstr << "Setting not found exception caught. Exchange striction unit cell lookup, setting llg.exch_striction.UnitCells" << i << " not set up properly for direction " << d << ". Setting " << snf.getPath() << " must be set.";
                                            std::string errstr=errsstr.str();
                                            error::errMessage(errstr);
                                        }
                                    }
                                    opsstr.str("");
                                    opsstr.clear();
                                    opsstr << "Unit cell direction " << uc;
                                    opstr=opsstr.str();

                                    FIXOUTVEC(config::Info,opstr,util::esuc(i,uc,0),util::esuc(i,uc,1),util::esuc(i,uc,2));
                                }
                                config::Info << std::endl;
                            }

                        }
                    }
                }
            }
            //for storing the field for each atom in the unit cell
            fields::sofield.resize(geom::ucm.GetNMS(),3);

            if(!config::cfg.exists("llg.staggered_field"))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Setting, llg.staggered)field not set, defaulting to zero.");
                for(unsigned int i = 0 ;  i < geom::ucm.GetNMS() ; i++)
                {
                    fields::sofield(i,0)=0.0;
                    fields::sofield(i,1)=0.0;
                    fields::sofield(i,2)=0.0;
                }
            }
            else
            {

                libconfig::Setting &sofset = config::cfg.lookup("llg.staggered_field");
                config::Info << "Staggered field Details:" << std::endl;
                config::Info << "========================" << std::endl;
                for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
                {
                    std::stringstream sofsstr;
                    sofsstr << "sapplied" << i+1;
                    std::string sofstr;
                    sofstr=sofsstr.str();
                    for(unsigned int j = 0 ; j < 3 ; j++)
                    {
                        try
                        {
                            fields::sofield(i,j)=sofset[sofstr.c_str()][j];
                        }
                        catch(const libconfig::SettingNotFoundException &snf)
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            std::stringstream errsstr;
                            errsstr << "Setting not found exception caught. Staggered field not set up properly. Setting " << snf.getPath() << " must be set.";

                            std::string errstr=errsstr.str();
                            error::errMessage(errstr);
                        }
                    }

                    FIXOUTVEC(config::Info,sofstr,fields::sofield(i,0),fields::sofield(i,1),fields::sofield(i,2));
                }

                config::Info << "========================" << std::endl;
            }
            config::Info << "Filling staggered field arrays..." << std::flush;
            for(unsigned int i = 0 ; i < geom::nspins ; i++)
            {

                unsigned int splu=geom::lu(i,3);
                fields::Hstagx(i)=fields::sofield(splu,0);
                fields::Hstagy(i)=fields::sofield(splu,1);
                fields::Hstagz(i)=fields::sofield(splu,2);
            }
            config::Info << "Done" << std::endl;

            if(!setting.lookupValue("dt",dt))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Time-step not set (llg.dt (double)), setting default 1e-16 [s]");
                dt=0.1e-15;
                FIXOUT(config::Info,"Timestep (default):" << dt << " seconds" << std::endl);
                rdt=dt*gyro;
                FIXOUT(config::Info,"Reduced timestep (default):" << rdt << std::endl);

            }
            else
            {
                FIXOUT(config::Info,"Timestep:" << dt << " seconds" << std::endl);
                rdt=dt*gyro;
                FIXOUT(config::Info,"Reduced timestep:" << rdt << std::endl);
            }
            if(!setting.lookupValue("intscheme",intschemestr))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Integration scheme not set (llg.intscheme (string)), setting default heun");
                intschemestr="heun";
                intscheme=0;
                FIXOUT(config::Info,"Numerical integration method (default):" << "Heun" << std::endl);
            }
            else
            {
                //Convert the string to lower case
                std::transform(intschemestr.begin(),intschemestr.end(),intschemestr.begin(),[](unsigned char c){return std::tolower(c);});
                //Check it is one of the recognised integration schemes
                if(intschemestr=="heun")
                {
                    intscheme=0;
                }
                else if(intschemestr=="rk4")
                {
                    intscheme=1;
                }
                else if(intschemestr=="rk5")
                {
                    intscheme=2;
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Read the integration scheme llg.intscheme (string) correctly, but the routine is not recognised.");
                }

                FIXOUT(config::Info,"Numerical integration method:" << intschemestr << std::endl);
            }
            int count=3;
            for(unsigned int i = 0 ; i < 3 ;i++)
            {
                try
                {
                    applied[i]=setting["applied"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errWarnPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " component " << i;
                    std::string errstr=errsstr.str();
                    error::errWarning(errstr);
                    count--;
                }
            }
            if(count==0)
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Field not specified, setting default to zero.");
                applied[0]=0.0;applied[1]=0.0;applied[2]=0.0;
                FIXOUTVEC(config::Info,"Applied field (default):",applied[0],applied[1],applied[2]);
            }
            else if(count==3)
            {
                 FIXOUTVEC(config::Info,"Applied field:",applied[0],applied[1],applied[2]);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Part of the field was not specified properly, check the warnings.");
            }

            if(!setting.lookupValue("update",spins::update))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read the spin update, setting default to 100 (timesteps)");
                spins::update=100;
                FIXOUT(config::Info,"Update of spins (default):" << spins::update << " [Timesteps]" << std::endl);
            }
            else
            {
                FIXOUT(config::Info,"Update of spins (default):" << spins::update << " [Timesteps]" << std::endl);
            }
            if(!setting.lookupValue("MagnetizationCalculationMethod:",spins::mag_calc_method))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read method for calculating magnetization (llg:MagnetizationCalculationMethod). Setting default to 0 (overall magnetization).");
                spins::mag_calc_method=0;
                FIXOUT(config::Info,"Magnetization calculatioin method (default):" << spins::mag_calc_method << std::endl);

            }
            else
            {
                FIXOUT(config::Info,"Magnetization calculatioin method:" << spins::mag_calc_method << std::endl);
            }
            if(spins::mag_calc_method==9 && geom::Nkset==true)
            {
                util::magDiscSize.resize(3);
                util::magDiscSize.IFill(0);
                for(unsigned int xyz = 0 ; xyz < 3 ; xyz++)
                {
                    try
                    {
                        util::magDiscSize[xyz]=setting["NoMagDisc"][xyz];
                    }
                    catch(const libconfig::SettingNotFoundException &snf)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        std::stringstream errsstr;
                        errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " component " << xyz;
                        std::string errstr=errsstr.str();
                        error::errMessage(errstr);
                    }
                }
                FIXOUTVEC(config::Info,"Number of intervals for calculating discrete magnetization:",util::magDiscSize[0],util::magDiscSize[1],util::magDiscSize[2]);

                double fractpart,intpart;

                double no_kpx_pd=static_cast<double>(geom::dim[0]*geom::Nk[0])/static_cast<double>(util::magDiscSize[0]);
                fractpart=std::modf(no_kpx_pd,&intpart);
                if(fabs(fractpart)>1e-12)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    printf ("No. points per discretization cell %f = %f + %f (the fractional part is bigger than 1e-12) \n", no_kpx_pd, intpart, fractpart);
                    error::errMessage("Error in discretization (in x). At the moment the code requires a perfect number of points so that geom::dim*geom::Nk/magDiscSize is EXACTLY an integer (to an accuracy of 1e-12).");
                }
                double no_kpy_pd=static_cast<double>(geom::dim[1]*geom::Nk[1])/static_cast<double>(util::magDiscSize[1]);
                fractpart=std::modf(no_kpy_pd,&intpart);
                if(fabs(fractpart)>1e-12)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    printf ("No. points per discretization cell %f = %f + %f (the fractional part is bigger than 1e-12) \n", no_kpy_pd, intpart, fractpart);
                    error::errMessage("Error in discretization (in y). At the moment the code requires a perfect number of points so that geom::dim*geom::Nk/magDiscSize is EXACTLY an integer (to an accuracy of 1e-12).");
                }
                double no_kpz_pd=static_cast<double>(geom::dim[2]*geom::Nk[2])/static_cast<double>(util::magDiscSize[2]);
                fractpart=std::modf(no_kpz_pd,&intpart);
                if(fabs(fractpart)>1e-12)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Error in discretization (in z). At the moment the code requires a perfect number of points so that geom::dim*geom::Nk/magDiscSize is EXACTLY an integer (to an accuracy of 1e-12).");
                    printf ("No. points per discretization cell %f = %f + %f (the fractional part is bigger than 1e-12) \n", no_kpz_pd, intpart, fractpart);
                }
                if(setting.lookupValue("DiscOutputFormat",util::disOutForm))
                {
                    FIXOUT(config::Info,"Output format for mag disc output (MagnetizationCalculationMethod):" << util::disOutForm << std::endl);
                }
                else
                {
                    error::errWarnPreamble(__FILE__,__LINE__);
                    error::errWarning("Could not read the DiscOutputFormat (llg:DiscOutputFormat (string)), defaulting to VTU.");
                    util::disOutForm="vtu";
                }
            }
            else if(spins::mag_calc_method==9 && geom::Nkset==false)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("To use magnetisation method 9 (discretization into cells) you have to define Nm.");
            }
            if(!setting.lookupValue("OutputMagnetization",spins::output_mag))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Output overall magnetization (llg.OutputMagnetization) not specified. Defaulting to true.");
            }
            FIXOUT(config::Info,"Initializing output of magnetization:" << std::endl);
            if(setting.lookupValue("SpinConfigMethod",scm))
            {
                FIXOUT(config::Info,"How are spins initially configured?:" << scm << std::endl);
                if(scm=="random")
                {
                    spins::setSpinsRandom();
                }
                else if(scm=="file")
                {
                    spins::setSpinsConfig();
                }
                else if(scm=="chequerx")
                {
                    spins::setSpinsChequerX();
                }
                else if(scm=="chequery")
                {
                    spins::setSpinsChequerY();
                }
                else if(scm=="chequerz")
                {
                    spins::setSpinsChequerZ();
                }
                else if(scm=="species")
                {
                    spins::setSpinsSpecies(setting);
                }
                else if(scm=="vampire")
                {
                    spins::setSpinsVampire(setting);
                }
                else if(scm=="restart_text")
                {
                    spins::setSpinsResText(setting);
                }
                else if(scm=="grid")
                {
                    spins::setSpinsGrid(setting);
                }
                else
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("Method for initialising spins not recognised (llg.SpinConfigMethod (string). Check you config file.");
                }
            }
            else
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read the method for applying the initial spin config. Setting llg.SpinConfigMethod (string), defaulting to file.");
                scm="file";
                spins::setSpinsConfig();
                FIXOUT(config::Info,"Initial spin configuratio (default):" << scm << std::endl);
            }
            if(!setting.lookupValue("OrdParTrans",trans))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read llg.OrdParTrans (int), setting default (0)");
                FIXOUT(config::Info,"Order parameter transform (default):" << trans << std::endl);
            }
            else
            {
                FIXOUT(config::Info,"Order parameter transform:" << trans << std::endl);
                if(trans==1)
                {
                    optrans.resize(geom::ucm.GetNMS(),3);
                    if(!config::cfg.exists("llg.OrdTrans"))
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("You requested to perform a transform to the order parameter calculation (e.g. to get the Neel vector) but you did not include the setting llg.OrdTrans in your config file.");
                    }
                    else
                    {
                        libconfig::Setting &subset = config::cfg.lookup("llg.OrdTrans");
                        for(unsigned int i = 0 ; i < geom::ucm.GetNMS() ; i++)
                        {
                            std::stringstream tsstr;
                            tsstr << "T" << i+1;
                            std::string tstr=tsstr.str();
                            for(unsigned int j = 0 ; j < 3 ; j++)
                            {
                                try
                                {
                                    optrans(i,j)=subset[tstr.c_str()][j];
                                }
                                catch(const libconfig::SettingNotFoundException &snf)
                                {
                                    std::stringstream errsstr;
                                    errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " component " << j << " for species " << i;
                                    std::string errstr=errsstr.str();
                                    error::errMessage(errstr);

                                }
                                if(fabs(optrans(i,j)-1.0>1e-6))
                                {
                                    error::errPreamble(__FILE__,__LINE__);
                                    error::errMessage("Each component of the spin transforms for outputting should have a magnitude of 1, check you input file under setting llg.OrdTrans");
                                }
                            }
                            std::stringstream opsstr;
                            opsstr << "Transformation for order parameter output for species " << i << ":";
                            std::string opstr=opsstr.str();
                            FIXOUTVEC(config::Info,opstr,optrans(i,0),optrans(i,1),optrans(i,2));

                        }
                    }
                }
            }
        }
        util::init_output();

        if(geom::ucm.NumAtomsUnitCell() > 5)
        {
            config::openLogFile();
        }

        config::printline(config::Info);
        //Output the prefactor to the LLG for each species to the output file
        for(unsigned int i = 0 ; i < geom::ucm.NumAtomsUnitCell() ; i++)
        {
            std::stringstream sstr,sstr1;
            std::string str=sstr.str();
            sstr << "Sigma prefactor for unit cell atom " << i << ":";
            sstr1 << "Prefactor for LLG for unit cell atom " << i << ":";
            geom::ucm.SetSigma(i,sqrt(2.0*kB*geom::ucm.GetDamping(i)/(geom::ucm.GetMu(i)*muB*dt*geom::ucm.GetGamma(i)*gyro)));
            geom::ucm.Setllgpf(i,-1./(1.0+geom::ucm.GetDamping(i)*geom::ucm.GetDamping(i)));

            if(i < 5)
            {
                FIXOUT(config::Info,str.c_str() << geom::ucm.GetSigma(i) << std::endl);
                str=sstr1.str();
                FIXOUT(config::Info,str.c_str() << geom::ucm.Getllgpf(i) << std::endl);
                config::printline(config::Info);
            }
            if(geom::ucm.NumAtomsUnitCell() > 5 && geom::logunit)
            {
                FIXOUT(config::Log,str.c_str() << geom::ucm.GetSigma(i) << std::endl);
                str=sstr1.str();
                FIXOUT(config::Log,str.c_str() << geom::ucm.Getllgpf(i) << std::endl);
                config::printline(config::Log);
            }
        }
        if(geom::ucm.NumAtomsUnitCell()>5)
        {
            for(unsigned int i = 0 ; i < 4 ; i++)
            {
                FIXOUT(config::Info,"             . . . " << "    . . ." << std::endl);
            }
            if(geom::logunit)
            {
                FIXOUT(config::Info,"FOR COMPLETE LLG INFO FOR UNIT CELL SEE LOG FILE:" << "   log.dat" << std::endl);
            }
            for(unsigned int i = 0 ; i < 4 ; i++)
            {
                FIXOUT(config::Info,"             . . . " << "    . . ." << std::endl);
            }

        }
        FIXOUT(config::Info,"Setting the llg and thermal prefactors to 1D arrays:" << std::flush);
        for(unsigned int i = 0 ; i < geom::nspins ; i++)
        {
            unsigned int aiuc=geom::lu(i,4);;
            geom::llgpf[i]=geom::ucm.Getllgpf(aiuc);
            geom::sigma[i]=geom::ucm.GetSigma(aiuc);
        }
        llg::Ts.resize(geom::ucm.GetNMS());
        llg::dps.resize(geom::ucm.GetNMS());
        llg::cps.resize(geom::ucm.GetNMS());

        SUCCESS(config::Info);

    }
    void integrate(unsigned int& t)
    {
        if(util::escheck && t%spins::update==0)
        {
            util::outESP << static_cast<double>(t)*dt;
            util::calc_es();
        }

#ifdef CUDA
        if(intscheme==0)
        {
            cullg::llgGPU(t);
        }
        else if(intscheme==1)
        {
            cullg::llgGPURK4(t);
        }
        else if(intscheme==2)
        {
            cullg::llgGPURK5(t);
        }
#else
        llgCPU::llgCPU(t);
#endif

    }
}
