// File: sim.h
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 23 May 2023 10:57:18 AM
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
    extern double rampinc;
}
#endif /*_SIM_H_*/
