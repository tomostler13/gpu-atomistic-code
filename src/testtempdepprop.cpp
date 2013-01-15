// File: geom.cpp
// Author:Tom Ostler
// Last-modified: 18 Dec 2012 18:03:53
//------header files for globals-------
#include "../inc/error.h"
#include "../inc/config.h"
#include "../inc/geom.h"
#include "../inc/mat.h"
#include "../inc/intmat.h"
#include "../inc/spins.h"
#include "../inc/fields.h"
#include "../inc/maths.h"
#include "../inc/util.h"
#include "../inc/neigh.h"
#include "../inc/LLBCPU.h"
#include "../inc/tdp.h"
#include "../inc/sim.h"
#include "../inc/nrf.h"
#include "../inc/mf.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#define FIXOUT(a,b) a.width(75);a << std::left << b;
void sim::testtempdepprop(int argc,char *argv[])
{
    assert(geom::gi);
    assert(tdp::tdpi);
    assert(fields::fi);
    assert(spins::si);
    assert(util::ui);
    assert(llb::LLBi);
    assert(mat::mi);
    assert(neigh::ni);
    assert(config::lcf);
    assert(sim::simi);
    assert(mf::mfi);

    if(tdp::toh!="uniform")
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Type of heating (toh) is not set to uniform, required for testing temperature dependent properties");
    }
    //min temperature
    double minTemp = 0.0000001;
    //maximum temperature
    double maxTemp = 700.0;
    //temperature step
    double dTemp = 5.0;
    //string containing output filename
    std::string opf;
    std::ofstream opfs;

    config::printline(config::Info);
    config::Info.width(45);config::Info << std::right << "*" << "**Simulation details***" << std::endl;
    FIXOUT(config::Info,"Simulation type:" << sim::sim_type << std::endl);

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

    libconfig::Setting &setting = config::cfg.lookup("sim.ttdp");
    try
    {
        setting.lookupValue("minTemp",minTemp);
        setting.lookupValue("maxTemp",maxTemp);
        setting.lookupValue("dTemp",dTemp);
        setting.lookupValue("outputFile",opf);
    }
    catch(const libconfig::SettingTypeException &stex)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Setting type error");
    }
    FIXOUT(config::Info,"Minimum temperature:" << minTemp << std::endl);
    FIXOUT(config::Info,"Maximum temperature:" << maxTemp << std::endl);
    FIXOUT(config::Info,"Change in temperature:" << dTemp << std::endl);
    FIXOUT(config::Info,"Output file:" << opf << std::endl);
    opfs.open(opf.c_str());
    if(!opfs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Cannot open file for writing");
    }
    opfs << "#T [K]\tChi Par [1/T]\tChi Perp [1/T]\tExchange Stiffness [J/m]\tEquilibrium magnetization [Red.]\n";
    for(double T = minTemp ; T < maxTemp ; T+=dTemp)
    {
        tdp::uniformtemp=T;
        tdp::Tcalc();
        tdp::mefp();
        tdp::chiparfp();
        tdp::chiperpfp();
        tdp::exchstifffp();
        opfs << tdp::systemp[0] << "\t" << tdp::syschipar[0] << "\t" << tdp::syschiperp[0] << "\t" << tdp::sysexchstiff[0] << "\t" << tdp::sysme[0] << std::endl;
    }
    opfs.close();
    if(opfs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errWarning("Could not close testing of temperature dependent properties file");
    }


}
