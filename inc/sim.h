// File: sim.h
// Author:Tom Ostler
// Created: 23 Jan 2013
// Last-modified: 23 Jan 2013 11:50:34
#include <iostream>
#include <fstream>
#ifndef _SIM_H_
#define _SIM_H_
namespace sim
{
    extern std::string sim_type;
    void initSim(int argc,char *argv[]);
	void MvT();
}
#endif /*_SIM_H_*/
