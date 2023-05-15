// File: sim.h
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 15 May 2023 03:45:54 PM
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
}
#endif /*_SIM_H_*/
