// File: sim.h
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 19 Aug 2022 13:13:54
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
}
#endif /*_SIM_H_*/
