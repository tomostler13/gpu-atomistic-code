// File: sim.h
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 31 May 2023 17:03:21
#include <iostream>
#include <fstream>
#ifndef _SIM_H_
#define _SIM_H_
namespace sim
{
    extern std::string sim_type;
    void initSim();
    void MvT();
    void suscep();
    void timeseries();
    void laser_heating();
    void ramp_field();
    void thermal_hyst();
    void simple();
    void ramp_stag();
    extern double rampinc,dHstg;
}
#endif /*_SIM_H_*/
