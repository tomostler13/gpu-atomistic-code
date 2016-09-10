// File: exch_routines.cpp
// Author: Tom Ostler
// Created: 05 Dec 2014
// Last-modified: 10 Sep 2016 17:00:32
// This source file was added to tidy up the file exch.cpp
// because it was becoming cumbersome to work with. The
// intention of this source file is to add a set of callable
// routines that do various neighbour list construction
// operations.
#include "../inc/arrays.h"
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/exch.h"
#include "../inc/intmat.h"
#include "../inc/llg.h"
#include "../inc/matrix_conv.h"
#include "../inc/defines.h"
#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>
namespace exch
{
    void readGlobalExch(int argc,char *argv[])
    {
        config::printline(config::Info);
        config::Info.width(45);config::Info << std::right << "*" << "**Exchange details***" << std::endl;
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
        libconfig::Setting &setting = config::cfg.lookup("exchange");
        if(setting.lookupValue("Read_exch_method",method))
        {
            FIXOUT(config::Info,"Method to read in exchange:" << method << std::endl);
        }
        else
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not read exchange method (exchange.Read_exchange_method (string)), check your config file.");
        }
        if(!setting.lookupValue("OutputExchange",outputJ))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Could not read exchange.OutputExchange (bool), defaulting to false.");
            outputJ=false;
        }
        if(!setting.lookupValue("OutputExchangeMatrix",oem))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Could not read whether to output the exchange matrix (exchange.OutputExchange (bool), defaulting to false.)");
            oem=false;
        }
        if(!setting.lookupValue("ReadExchangeMatrix",rem))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Could not determine whether exchange matrix should be read (exchange.ReadExchangeMatrix (bool)), defaulting to false.");
            rem=false;
        }
        if(!setting.lookupValue("ReadExchangeMatrix4Spin",rem4s))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Could not read whether you want to read the four spin CSR matrix. Defaulting to false.");
            rem4s=false;
        }
        else
        {
            FIXOUT(config::Info,"Reading 4 spin exchange matrix?:" << config::isTF(rem4s) << std::endl);
        }

        if(!setting.lookupValue("ExitAfterExchangeMatrix",eaem))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Could not read whether you want to exit after outputting the exchange matrix. exchange.ExitAfterExchangeMatrix (bool), defaulting to false.");
            eaem=false;
        }
        FIXOUT(config::Info,"Exit after outputting exchange matrix?:" << config::isTF(eaem) << std::endl);
        if(oem && rem)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("There is no point in reading and writing the exchange matrix because if you are reading it in you already have it and it would remain unchanged.");
        }
        FIXOUT(config::Info,"Output the exchange matrix in the current format:" << config::isTF(oem) << std::endl);
        FIXOUT(config::Info,"Read the exchange matrix:" << config::isTF(rem) << std::endl);
        if(rem)
        {
            if(setting.lookupValue("ExchangeMatrixFilename",exchMatFN))
            {
                FIXOUT(config::Info,"Exchange matrix filename:" << exchMatFN << std::endl);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Error reading filename for csr exchange matrix. Check setting exchange.csrMatrixName");
            }
        }
        FIXOUT(config::Info,"Output of exchange matrix:" << config::isTF(outputJ) << std::endl);
        FIXOUT(config::Info,"Reading exchange interactions from:" << "external file" << std::endl);
        if(!setting.lookupValue("Exchfile",readFile))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("You must specify an exchange file.");
        }
        FIXOUT(config::Info,"File:" << readFile << std::endl);
        if(method=="permute" || method=="mapint" || method=="unitcell")
        {
            try
            {
                exchcfg.readFile(readFile.c_str());
            }
            catch(const libconfig::FileIOException &fioex)
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("I/O error while reading exchange config file");
            }
            catch(const libconfig::ParseException &pex)
            {
                error::errPreamble(__FILE__,__LINE__);
                std::cerr << ". Parse error at " << pex.getFile()  << ":" << pex.getLine() << "-" << pex.getError() << "***\n" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        libconfig::Setting &GlobExch = exchcfg.lookup("exchange");
        if(!GlobExch.lookupValue("FourSpin",inc4spin))
        {
            error::errWarnPreamble(__FILE__,__LINE__);
            error::errWarning("Could not read if you want to include the four spin term (exchange.FourSpin (bool), defaulting to false.");
        }
        FIXOUT(config::Info,"Include four spin term?:" << config::isTF(inc4spin) << std::endl);
        if(method=="permute")
        {
            if(!GlobExch.lookupValue("MaxShells",max_shells))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not read the max number of shells of exchange to read in (exchange.MaxShells (int)). Must be set.");
            }
            evs.resize(3);
            evs[0]=GlobExch["Scale"][0];
            evs[1]=GlobExch["Scale"][1];
            evs[2]=GlobExch["Scale"][2];
            FIXOUTVEC(config::Info,"Scaling factor for exchange vectors:",evs[0],evs[1],evs[2]);
            numint.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells);
            exchvec.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_shells,3);
            shell_list.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS());
            numint.IFill(0);exchvec.IFill(0);J.IFill(0);shell_list.IFill(0);
        }
        else if(method=="mapint")//if this is true we no longer care about shells, the interactions are hard coded (based on the mesh points)
        {
            if(!GlobExch.lookupValue("MaxInteractions",max_int))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not read the max number of interactions of exchange to read in (exchange.MaxInteractions (int)).");
            }
            if(!GlobExch.lookupValue("TruncateExchange",cutexch))
            {
                error::errWarnPreamble(__FILE__,__LINE__);
                error::errWarning("Could not read if you wanted to trucate the interaction range (exchange.TruncateExchange (bool)). Defaulting to false");
                cutexch=false;
            }
            evs.resize(3);
            for(unsigned int i = 0 ; i < 3 ; i++)
            {
                try
                {
                    evs[i]=GlobExch["Scale"][i];
                }
                catch(const libconfig::SettingNotFoundException &snf)
                {
                    error::errPreamble(__FILE__,__LINE__);
                    std::stringstream errsstr;
                    errsstr << "Setting not found exception caught. Setting " << snf.getPath() << " must be set for mapint method.";
                    std::string errstr=errsstr.str();
                    error::errMessage(errstr);
                }
            }
            FIXOUTVEC(config::Info,"Scaling factor for exchange vectors:",evs[0],evs[1],evs[2]);
            FIXOUT(config::Info,"Truncate exchange?:" << config::isTF(cutexch) << std::endl);
            if(cutexch==true)
            {
                if(!GlobExch.lookupValue("TruncateRCut",rcut))
                {
                    error::errPreamble(__FILE__,__LINE__);
                    error::errMessage("You stated that you wanted to cut the interaction range but have not specified a cut-off radius (exchange.TruncateRCut (double) in metres).");
                }
                else
                {
                    FIXOUT(config::Info,"Cut-off radius:" << rcut << " [m]" << std::endl);
                }
            }
            FIXOUT(config::Info,"Maximum number of interactions:" << max_int << std::endl);
            numint.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_int);
            exchvec.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_int,3);
            shell_list.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS());
            numint.IFill(0);exchvec.IFill(0);shell_list.IFill(0);
        }
        else if(method=="unitcell")
        {
            if(!GlobExch.lookupValue("MaxInteractions",max_int))
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not read exchange.MaxInteractions (int). Must be set.");
            }
            FIXOUT(config::Info,"Maximum number of interactions:" << max_int << std::endl);
            numint.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_int);
            exchvec.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS(),max_int,3);
            shell_list.resize(geom::ucm.GetNMS(),geom::ucm.GetNMS());
            numint.IFill(0);exchvec.IFill(0);shell_list.IFill(0);
        }
        if(inc4spin)
        {
            if(GlobExch.lookupValue("MaxQuartetPermutations",max_4s))
            {
                FIXOUT(config::Info,"Maximum number of quartet permutations:" << max_4s << std::endl);
            }
            else
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("Could not read the setting exchange.MaxQuartetPermutations (int)");
            }
            read4spin();
        }
    }

    void permute(int argc,char *argv[])
    {
        get_exch_permute(argc,argv);
        checkdiag.resize((geom::nspins*2)-1);
        checkdiag.IFill(0);
        //then add the exchange constants to the matrix or to the interaction matrix
        if(config::exchm>0)//We calculate the exchange via a matrix multiplication (or part of it at least)
        {
            if(config::exchm==2)
            {
                exchpermute();
                //output the sparse matrix
                if(oem)
                {
                    if(config::offdiag==false)
                    {

                        std::ofstream opem("csr_exch_mat.dat");
                        if(!opem.is_open())
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("Could not open file for outputting dia exchange matrix.");
                        }
                        opem << geom::nspins << std::endl;
                        opem << adjncy.size() << std::endl;
                        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            opem << xadj[i] << "\t";
                        }
                        opem << xadj[geom::nspins] << std::endl;
                        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
                            {
                                opem << std::setprecision(16) << adjncy[j] << "\t" << dataxx[j] << "\t" << datayy[j] << "\t" << datazz[j] << "\t";
                            }
                        }
                        opem.close();
                        if(opem.is_open())
                        {
                            error::errWarnPreamble(__FILE__,__LINE__);
                            error::errWarning("Could not close the file for outputting the CSR exchange matrix.");
                        }

                    }
                    else
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("FUTURE IMPLEMENTATION: At the moment the outputting of the off-diagonal CSR format is not possible.");
                    }

                }
                //For debugging the CSR neighbour list
                /*for(unsigned int i = 0 ; i < geom::nspins ; i++)
                  {
                  std::cout << "Atom " << i << " has neighbours" << std::endl;
                  for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
                  {
                  std::cout << adjncy[j] << "\t" << dataxx[j] << "\t" << std::flush;
                  }
                  std::cout << std::endl;
                  }
                  std::cin.get();*/

                if(config::exchm==1)
                {//then convert CSR to DIA

                    //call the routine to convert the CSR to a DIA sparse matrix format
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_diag(geom::nspins,xadj,adjncy,diagoffset,diagnumdiag,checkdiag,dataxx,datayy,datazz);
                    SUCCESS(config::Info);
                    FIXOUT(config::Info,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                    //for debugging the data array
                    /*for(int i = 0 ; i < diagnumdiag ; i++)
                      {
                      int os=diagoffset[i];
                      std::cout << "Offset " << i << " is " << os << std::endl;
                      for(int j = 0 ; j < geom::nspins ; j++)
                      {
                      std::cout << j << "\t" << i*geom::nspins+j << "\t" << datazz[i*geom::nspins+j] << std::endl;
                      }
                      std::cin.get();
                      }*/
                    if(oem)
                    {
                        if(config::offdiag)
                        {
                            std::ofstream opem("dia_exch_mat.dat");
                            if(!opem.is_open())
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("Could not open file for outputting dia exchange matrix.");
                            }
                            opem << geom::nspins << std::endl;
                            opem << diagnumdiag << std::endl;
                            for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                            {
                                opem << diagoffset[i] << "\t";
                            }
                            opem << diagoffset[diagnumdiag-1];
                            for(unsigned int i = 0 ; i < diagnumdiag ; i++)
                            {
                                for(unsigned int j = 0 ; j < geom::nspins ; j++)
                                {
                                    opem << dataxx[i*geom::nspins+j] << "\t" << datayy[i*geom::nspins+j] << "\t" << datazz[i*geom::nspins+j] << "\t";
                                }
                            }
                            opem.close();
                            if(opem.is_open())
                            {
                                error::errWarnPreamble(__FILE__,__LINE__);
                                error::errWarning("Could not close the file for outputting the DIA exchange matrix.");
                            }
                        }
                        else
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("FUTURE IMPLEMENTATION: At the moment outputting of the off-diagonal in DIA format is not possible.");
                        }

                    }

                }
#ifdef CUDA
                if(config::exchm>98)
                {
                    //call the routine to convert the CSR to a DIA sparse matrix format
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_diag(geom::nspins,xadj,adjncy,diagoffset,diagnumdiag,checkdiag,dataxx,datayy,datazz);
                    SUCCESS(config::Info);
                    FIXOUT(config::Info,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                }
                if(config::offdiag && config::exchm>98)
                {
                    offdiagoffset.resize(diagoffset.size());
                    for(unsigned int i = 0 ; i < offdiagoffset.size() ; i++)
                    {
                        offdiagoffset[i]=diagoffset[i];
                    }
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_offdiag(geom::nspins,xadj,adjncy,offdiagoffset,diagnumdiag,checkdiag,dataxy,dataxz,datayx,datayz,datazx,datazy);
                    SUCCESS(config::Info);
                    FIXOUT(config::Log,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                    FIXOUT(config::Log,"Total number of non-zero diagonals for the anti-symmetric part of the exchange tensor:" << offdiagnumdiag << std::endl);
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for off-diagonal (antisym) part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < offdiagnumdiag-1 ; i++)
                    {
                        config::Log << offdiagoffset[i] << ",";
                    }
                    config::Log << offdiagoffset[offdiagnumdiag-1] << " ] " << std::endl;
                    if(oem)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("FUTURE IMPLEMENTATION: At the moment it is not possible to output the exchange matrix of the hybrid format.");
                    }


                }
#endif /*CUDA*/
                if(config::offdiag && config::exchm==1)
                {
                    offdiagoffset.resize(diagoffset.size());
                    for(unsigned int i = 0 ; i < offdiagoffset.size() ; i++)
                    {
                        offdiagoffset[i]=diagoffset[i];
                    }
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_offdiag(geom::nspins,xadj,adjncy,offdiagoffset,diagnumdiag,checkdiag,dataxy,dataxz,datayx,datayz,datazx,datazy);
                    SUCCESS(config::Info);
                    FIXOUT(config::Log,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                    FIXOUT(config::Log,"Total number of non-zero diagonals for the anti-symmetric part of the exchange tensor:" << offdiagnumdiag << std::endl);
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for off-diagonal (antisym) part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < offdiagnumdiag-1 ; i++)
                    {
                        config::Log << offdiagoffset[i] << ",";
                    }
                    config::Log << offdiagoffset[offdiagnumdiag-1] << " ] " << std::endl;


                }
                if(config::exchm!=2)//then we are not using CSR
                {
                    xadj.clear();
                    adjncy.clear();
                    checkdiag.clear();
                }
            }
        }
        else if(config::exchm==0)//add the exchange to the interaction matrix
        {
            fft();
        }


        if(config::exchm>98)
        {
            hybrid();
        }//end of else if(config::exchm>98)
    }

    void mapint(int argc,char *argv[])
    {
        get_exch_mapint(argc,argv);
        checkdiag.resize((geom::nspins*2)-1);
        checkdiag.IFill(0);
        if(config::exchm>0)//We calculate the exchange via a matrix multiplication (or part of it at least)
        {

            if(config::exchm==2)
            {
                exchmapint();
                //output the spare matrix
                if(oem)
                {
                    if(config::offdiag==false)
                    {

                        std::ofstream opem("csr_exch_mat.dat");
                        if(!opem.is_open())
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("Could not open file for outputting dia exchange matrix.");
                        }
                        opem << geom::nspins << std::endl;
                        opem << adjncy.size() << std::endl;
                        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            opem << xadj[i] << "\t";
                        }
                        opem << xadj[geom::nspins] << std::endl;
                        for(unsigned int i = 0 ; i < geom::nspins ; i++)
                        {
                            for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
                            {
                                opem << adjncy[j] << "\t" << dataxx[j] << "\t" << datayy[j] << "\t" << datazz[j] << "\t";
                            }
                        }
                        opem.close();
                        if(opem.is_open())
                        {
                            error::errWarnPreamble(__FILE__,__LINE__);
                            error::errWarning("Could not close the file for outputting the CSR exchange matrix.");
                        }

                    }
                    else
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("FUTURE IMPLEMENTATION: At the moment the outputting of the off-diagonal CSR format is not possible.");
                    }

                }
                //For debugging the CSR neighbour list
                /*for(unsigned int i = 0 ; i < geom::nspins ; i++)
                  {
                  std::cout << "Atom " << i << " has neighbours" << std::endl;
                  for(unsigned int j = xadj[i] ; j < xadj[i+1] ; j++)
                  {
                  std::cout << adjncy[j] << "\t" << dataxx[j] << "\t" << std::flush;
                  }
                  std::cout << std::endl;
                  }
                  std::cin.get();*/

                if(config::exchm==1)
                {//then convert CSR to DIA

                    //call the routine to convert the CSR to a DIA sparse matrix format
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_diag(geom::nspins,xadj,adjncy,diagoffset,diagnumdiag,checkdiag,dataxx,datayy,datazz);
                    SUCCESS(config::Info);
                    FIXOUT(config::Info,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                    //for debugging the data array
                    /*for(int i = 0 ; i < diagnumdiag ; i++)
                      {
                      int os=diagoffset[i];
                      std::cout << "Offset " << i << " is " << os << std::endl;
                      for(int j = 0 ; j < geom::nspins ; j++)
                      {
                      std::cout << j << "\t" << i*geom::nspins+j << "\t" << datazz[i*geom::nspins+j] << std::endl;
                      }
                      std::cin.get();
                      }*/
                    if(oem)
                    {
                        if(config::offdiag)
                        {
                            std::ofstream opem("dia_exch_mat.dat");
                            if(!opem.is_open())
                            {
                                error::errPreamble(__FILE__,__LINE__);
                                error::errMessage("Could not open file for outputting dia exchange matrix.");
                            }
                            opem << geom::nspins << std::endl;
                            opem << diagnumdiag << std::endl;
                            for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                            {
                                opem << diagoffset[i] << "\t";
                            }
                            opem << diagoffset[diagnumdiag-1];
                            for(unsigned int i = 0 ; i < diagnumdiag ; i++)
                            {
                                for(unsigned int j = 0 ; j < geom::nspins ; j++)
                                {
                                    opem << dataxx[i*geom::nspins+j] << "\t" << datayy[i*geom::nspins+j] << "\t" << datazz[i*geom::nspins+j] << "\t";
                                }
                            }
                            opem.close();
                            if(opem.is_open())
                            {
                                error::errWarnPreamble(__FILE__,__LINE__);
                                error::errWarning("Could not close the file for outputting the DIA exchange matrix.");
                            }
                        }
                        else
                        {
                            error::errPreamble(__FILE__,__LINE__);
                            error::errMessage("FUTURE IMPLEMENTATION: At the moment outputting of the off-diagonal in DIA format is not possible.");
                        }

                    }

                }
#ifdef CUDA
                if(config::exchm>98)
                {
                    //call the routine to convert the CSR to a DIA sparse matrix format
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_diag(geom::nspins,xadj,adjncy,diagoffset,diagnumdiag,checkdiag,dataxx,datayy,datazz);
                    SUCCESS(config::Info);
                    FIXOUT(config::Info,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                }
                if(config::offdiag && config::exchm>98)
                {
                    offdiagoffset.resize(diagoffset.size());
                    for(unsigned int i = 0 ; i < offdiagoffset.size() ; i++)
                    {
                        offdiagoffset[i]=diagoffset[i];
                    }
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_offdiag(geom::nspins,xadj,adjncy,offdiagoffset,diagnumdiag,checkdiag,dataxy,dataxz,datayx,datayz,datazx,datazy);
                    SUCCESS(config::Info);
                    FIXOUT(config::Log,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                    FIXOUT(config::Log,"Total number of non-zero diagonals for the anti-symmetric part of the exchange tensor:" << offdiagnumdiag << std::endl);
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for off-diagonal (antisym) part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < offdiagnumdiag-1 ; i++)
                    {
                        config::Log << offdiagoffset[i] << ",";
                    }
                    config::Log << offdiagoffset[offdiagnumdiag-1] << " ] " << std::endl;
                    if(oem)
                    {
                        error::errPreamble(__FILE__,__LINE__);
                        error::errMessage("FUTURE IMPLEMENTATION: At the moment it is not possible to output the exchange matrix of the hybrid format.");
                    }


                }
#endif /*CUDA*/
                if(config::offdiag && config::exchm==1)
                {
                    offdiagoffset.resize(diagoffset.size());
                    for(unsigned int i = 0 ; i < offdiagoffset.size() ; i++)
                    {
                        offdiagoffset[i]=diagoffset[i];
                    }
                    FIXOUT(config::Info,"Converting CSR to DIA (exchange tensor diagonals):" << std::flush);
                    matconv::csr_to_dia_offdiag(geom::nspins,xadj,adjncy,offdiagoffset,diagnumdiag,checkdiag,dataxy,dataxz,datayx,datayz,datazx,datazy);
                    SUCCESS(config::Info);
                    FIXOUT(config::Log,"Total number of non-zero diagonals (size of offset array):" << diagnumdiag << std::endl);
                    config::openLogFile();
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for diagonal part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < diagnumdiag-1 ; i++)
                    {
                        config::Log << diagoffset[i] << ",";
                    }
                    config::Log << diagoffset[diagnumdiag-1] << " ] " << std::endl;
                    FIXOUT(config::Log,"Total number of non-zero diagonals for the anti-symmetric part of the exchange tensor:" << offdiagnumdiag << std::endl);
                    config::printline(config::Log);
                    FIXOUT(config::Log,"Outputting offset for off-diagonal (antisym) part of interaction matrix:" << std::endl);
                    config::Log << " [ ";
                    for(unsigned int i = 0 ; i < offdiagnumdiag-1 ; i++)
                    {
                        config::Log << offdiagoffset[i] << ",";
                    }
                    config::Log << offdiagoffset[offdiagnumdiag-1] << " ] " << std::endl;


                }
                if(config::exchm!=2)//then we are not using CSR
                {
                    xadj.clear();
                    adjncy.clear();
                    checkdiag.clear();
                }
            }
        }
        else if(config::exchm==0)//add the exchange to the interaction matrix
        {
            directfft();
        }

        if(config::exchm>98)
        {
            directhybrid();
        }//end of else if(config::exchm>98)
    }
}
