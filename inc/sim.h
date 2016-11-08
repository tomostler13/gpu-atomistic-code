// File: sim.h
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 18 Oct 2016 13:14:48
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
}
#endif /*_SIM_H_*/
