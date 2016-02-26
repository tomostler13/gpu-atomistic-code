#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "../../../inc/array.h"
#include "../../../inc/array2d.h"
#include "../../../inc/array3d.h"
#include "../../../inc/array4d.h"
#include "../../../inc/array5d.h"
#include "../../../inc/defines.h"
#include "../../../inc/error.h"
int isInteger(double x)
{
    if(std::ceil(x)!=std::floor(x))
    {
        return(1);
    }
    else
    {
        return(0);
    }
}
int main(int argc,char *argv[])
{
    if(argc<8)
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Incorrect usage, should be: ./convSPRKKR <inputfile> <output exchange file name> <scale x> <scale y> <scale z> <scale Jij> <units (eV or Ry)>");
    }
    std::string fnm=argv[1],opn=argv[2],units=argv[7];
    double ns[3]={atof(argv[3]),atof(argv[4]),atof(argv[5])};
    double scaleJ = atof(argv[6]);
    std::ofstream logf("log.dat"),opf(opn.c_str());
    if(!logf.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open log.dat file.");
    }
    if(!opf.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open output exchange file, check the filename.");
    }
    std::ifstream ifs(fnm.c_str());
    if(!ifs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not open input file.");
    }
    std::string sdump,inst;
    //Get and dump the header
    ifs >> sdump >> sdump >> sdump >> sdump >> sdump >> sdump;
    ifs >> sdump >> sdump >> sdump >> sdump >> sdump;
    std::getline(ifs,inst);
    FIXOUT(logf,"Opened input file:" << fnm << std::endl);
    FIXOUT(logf,"Opened file for outputting exchange:" << "Continuing" << std::endl);
    FIXOUTVEC(logf,"Scaling vector for lookup coordinates:",ns[0],ns[1],ns[2]);
    FIXOUT(logf,"Scaling factor for Jij values:" << scaleJ << std::endl);
    if(units=="eV" || units=="Ry")
    {
        FIXOUT(logf,"Units of exchange:" << units << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The units were not recognised. Either choose eV or Ry (take note of upper/lower case, it matters here.");
    }
    bool brk=false;
    int nl=0;
    std::string search="IQ";
    while(brk==false)
    {
        std::getline(ifs,inst);
    //    std::cout << inst << std::endl;
        nl++;
        if(inst.find(search)!=std::string::npos)
        {
            brk=true;
        }
        if(nl>1e6)
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The number of lines read before reading of the exchange exceeded 1,000,000. This is probably an error.");
        }
    }
    unsigned numsite=nl-2;
    FIXOUT(logf,"Number of magnetic sites:" << numsite << std::endl);
    //should not have reached the end of the file so don't (in theory) have to clear to eof flag
    //with ifile.clear()
    ifs.seekg(0,std::ios::beg);
    //Get and dump the header
    ifs >> sdump >> sdump >> sdump >> sdump >> sdump >> sdump;
    ifs >> sdump >> sdump >> sdump >> sdump >> sdump;
    std::getline(ifs,inst);
    //store the number of atoms on that site (for CPA) - at the moment this must be the same for each site
    Array<int> naps;
    naps.resize(numsite);naps.IFill(0);
    //store the id's of the atoms on the site (naps)
    Array2D<int> aid;
    //store the compositions
    Array2D<double> comp;
    //store the numbers
    Array2D<int> numbers;
    for(unsigned int i = 0 ; i < numsite ; i++)
    {
        int site=0;
        ifs >> site;
        if(site!=(i+1))
        {
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("The number of atoms is not correct.");
        }
        ifs >> naps[i];
        logf << "Number of atoms at site " << i+1 << " = " << naps[i] << std::endl;
        if(i>0)
        {
            if(naps[i]!=naps[i-1])
            {
                error::errPreamble(__FILE__,__LINE__);
                error::errMessage("The number of atoms per site MUST be the same for each atomic site at the moment. It has been detected that it is not. Please check your SPRKKR file.");
                //Note - if you want to have a different number of atoms per site (not sure what this means for ASD) then you will probably want to replace
                //the Array2D (aid) with a std::vector< std::vector<int> > that can have arbitary dimensions on each of the first dimensions.
            }
        }
        else
        {
            aid.resize(numsite,naps[0]);
            aid.IFill(0);
            comp.resize(numsite,naps[0]);
            comp.IFill(0);
            numbers.resize(numsite,naps[0]);
            numbers.IFill(0);
        }
        //for debuggin
        //std::cout << i+1 << "\t" << naps[i];
        //read the number and it's associated weight
        for(unsigned int j = 0 ; j < naps[i] ; j++)
        {
            ifs >> numbers(i,j) >> comp(i,j);
            //for debugging
            logf << "\t" << numbers(i,j) << "\t" << comp(i,j) << std::endl;
        }
        //for debugging
        //std::cout << std::endl;
    }
    int totnoint=0;
    bool check=false;
    //we don't know how big the interaction range is yet, so pre-read it
    while(check==false)
    {
        int IQ=-1,IT=-1,JQ=-1,JT=-1;
        ifs >> inst >> sdump >> IQ >> sdump >> sdump >> IT >> sdump >> sdump >> JQ >> sdump >> sdump >> JT;

        if(IQ>numsite/2)
        {
            check=true;
        }

        if(check==false)
        {
            //don't care about the next line
            std::getline(ifs,sdump);//std::getline(ifs,sdump);
            //debug
            //std::cout << IQ << "\t" << IT << "\t" << JQ << "\t" << JT << std::endl;
            std::getline(ifs,sdump);std::getline(ifs,sdump);
            double x,y,z,dist,JRy,JeV;
            ifs >> sdump >> sdump >> sdump >> sdump >> sdump >> x >> y >> z >> dist >> JRy >> JeV;
            //debug
            //std::cout << x << "\t" << y << "\t" << z << "\t" << dist << "\t" << JRy << "\t" << JeV << std::endl;//std::cin.get();


            totnoint++;
        }
    }
    FIXOUT(logf,"Total number of interactions:" << totnoint << std::endl);
    FIXOUT(logf,"Maximum number of interactions between any two species:" << totnoint/(numsite*naps[0]) << std::endl);
    //Maximum interactions should be an integer (i.e. the same number for each species)
    double tmi=static_cast<double>(totnoint)/static_cast<double>(numsite*naps[0]);
    if(isInteger(tmi)==0)
    {
        FIXOUT(logf,"Number of interactions between any two species was found to be an integer:" << "Continuing" << std::endl);
    }
    else
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("The number of interactions between any species was not an integer");
    }

    unsigned int mi=static_cast<double>(tmi);
    Array<double> JijRy,JijeV;
    Array2D<int> speclu;
    Array2D<double> vec;

    JijeV.resize(totnoint);
    JijeV.IFill(0.0);
    JijRy.resize(totnoint);
    JijRy.IFill(0.0);
    vec.resize(totnoint,3);
    vec.IFill(0.0);
    //element 0 is the id of speci and 1 is the id of specj
    speclu.resize(totnoint,2);
    speclu.IFill(-1);

    //no read the file again
    ifs.seekg(0,std::ios::beg);
    //Get and dump the header
    ifs >> sdump >> sdump >> sdump >> sdump >> sdump >> sdump;
    ifs >> sdump >> sdump >> sdump >> sdump >> sdump;
    std::getline(ifs,inst);
    for(unsigned int i = 0 ; i < numsite ; i++)
    {
        int site=0;
        ifs >> sdump;
        ifs >> sdump;
        for(unsigned int j = 0 ; j < naps[i] ; j++)
        {
            ifs >> sdump >> sdump;
        }
    }
    //Loop over all of the interactions
    for(unsigned int i = 0 ; i < totnoint ; i++)
    {
        int speci=-1,specj=-1;
        int IQ=-1,IT=-1,JQ=-1,JT=-1;
        ifs >> inst >> sdump >> IQ >> sdump >> sdump >> IT >> sdump >> sdump >> JQ >> sdump >> sdump >> JT;
        //species of i
        speci=IT%naps[0];
        //species of j
        specj=JT%naps[0];
        std::getline(ifs,sdump);
        std::getline(ifs,sdump);std::getline(ifs,sdump);
        double x,y,z,dist,JRy,JeV;
        ifs >> sdump >> sdump >> sdump >> sdump >> sdump >> x >> y >> z >> dist >> JRy >> JeV;
        vec(i,0)=x;
        vec(i,1)=y;
        vec(i,2)=z;

        speclu(i,0)=speci;
        speclu(i,1)=specj;
        if(speclu(i,0)==-1)
        {
            std::cerr << IT << "\t" << naps[0] << "\t" << JT%naps[0] << std::endl;
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not get species (i) id.");
        }
        if(speclu(i,1)==-1)
        {
            std::cerr << JT << "\t" << naps[0] << "\t" << JT%naps[0] << std::endl;
            error::errPreamble(__FILE__,__LINE__);
            error::errMessage("Could not get species (j) id.");
        }
        JijeV(i)=JeV;
        JijRy(i)=JRy;

    }
    opf << "exchange:\n{\n\tFourSpin = FALSE;\n\tScale = [ 1.0 , 1.0 , 1.0 ];" << std::endl;
    opf << "\tMaxInteractions = " << mi << ";" << std::endl;
    opf << "\tTruncateExchange = FALSE;\n};" << std::endl;
    for(unsigned int i = 0 ; i < naps[0] ; i++)
    {
        for(unsigned int j = 0 ; j < naps[0] ; j++)
        {

            int loccount=0;
            std::stringstream sstr;
            sstr << "exchange_" << i << "_" << j << ":";
            std::string str=sstr.str();
            opf << str << std::endl;
            opf << "{\n\tNum_Interactions = " << mi << ";\n\tunits = \"" << units << "\";\n";
            std::stringstream opfsstr;
            opfsstr << "rJ_" << i << "_" << j << ".dat";
            std::string opfstr=opfsstr.str();
            std::ofstream raw(opfstr.c_str());
            //now output the interaction between each species
            //loop over all of the read interactions and check if they are between each species
            for(unsigned int l = 0 ; l < totnoint/2 ; l++)
            {
//            std::cout << "l= " << l << "\tlocout = " << loccount << std::endl;
                //if(l<mi)
                {
                    if(speclu(l,0)==i &&  speclu(l,1)==j)
                    {
                        opf << std::setprecision(10) << std::fixed << "\tInteraction" << loccount+1 << "Vec = [ " << vec(l,0)*ns[0] << " , " << vec(l,1)*ns[1] << " , " << vec(l,2)*ns[2] << " ];\n";

                        if(units=="eV")
                        {
                            opf << "\tJ" << loccount+1 << "_1 = [ " << JijeV(l)*scaleJ << " , 0.0000 , 0.0000 ];\n";
                            opf << "\tJ" << loccount+1 << "_2 = [ 0.0000 , " << JijeV(l)*scaleJ << " , 0.0000 ];\n";
                            opf << "\tJ" << loccount+1 << "_3 = [ 0.0000 , 0.0000 , " << JijeV(l)*scaleJ << " ];\n";
                        }
                        else if(units=="Ry")
                        {
                            opf << "\tJ" << loccount+1 << "_1 = [ " << JijRy(l)*scaleJ << " , 0.0000 , 0.0000 ];\n";
                            opf << "\tJ" << loccount+1 << "_2 = [ 0.0000 , " << JijRy(l)*scaleJ << " , 0.0000 ];\n";
                            opf << "\tJ" << loccount+1 << "_3 = [ 0.0000 , 0.0000 , " << JijRy(l)*scaleJ << " ];\n";
                        }
                        raw << vec(l,0)*ns[0] << "\t" << vec(l,1)*ns[1] << "\t" << vec(l,2)*ns[2] << "\t" << JijeV(l)*scaleJ << std::endl;
                        loccount++;

                    }
                }
            }
            raw.close();
            if(loccount!=mi)
            {
                error::errPreamble(__FILE__,__LINE__);
//                std::cout <<"\n"<< loccount << "\t" << mi << "\t" << totnoint << std::endl;
                error::errMessage("Inconsistency in the number of interactions found");
            }
            opf << "};\n";
        }
    }
    ifs.close();
    if(ifs.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not close the SPRKKR (input) file.");
    }
    logf.close();
    if(logf.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not close the log.dat (output) file.");
    }
    opf.close();
    if(opf.is_open())
    {
        error::errPreamble(__FILE__,__LINE__);
        error::errMessage("Could not close the exchange output file. Check it's contents.");
    }
    return(0);
}
