// File: sim.h
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 24 Nov 2014 14:04:49
#include <iostream>
#include <fstream>
#ifndef _SIM_H_
#define _SIM_H_
namespace sim
{
    extern std::string sim_type;
    void initSim(int argc,char *argv[]);
    void MvT(int argc,char *argv[]);
    void suscep(int argc,char *argv[]);
    void timeseries(int argc,char *argv[]);
    void laser_heating(int argc,char *argv[]);
}
#endif /*_SIM_H_*/
