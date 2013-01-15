#include <iostream>
#include <fstream>
#ifndef _SIM_H_
#define _SIM_H_
#define FIXOUT(a,b) a.width(75);a << std::left << b;
namespace sim
{
    extern bool simi;
    extern double dt;
    extern double rdt;
    extern double llbpf;
    extern std::string sim_type;
    void initSim(int argc,char *argv[]);
    void testtempdepprop(int argc,char *argv[]);
    void testsuite(int argc,char *argv[]);
}
#endif /*_SIM_H_*/
