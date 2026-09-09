#include "common.hpp"
#include "init.hpp"
#include "init_interp.hpp"
#include "params.hpp"
#include "sigma_LO.hpp"
#include "sigma_NLO.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace params;

//extern "C"
//{
    // ih: input hadron, 1=pion, 2=kaon, 3=proton, 4=charged hadrons
    // ic: charge, 0=0, 1=+, -1=-
    // parton:
        //0    1    2    3    4    5    6    7    8     9    10
        //g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
    // io: order, 0=LO, 1=NLO
    // NOTE: cbar and bbar not in DSS, returns c or b instead
    // result is the fragmentation function (X factor divided out)

    /*void fdss_ (int &ih, int &ic, int &io, double &x, double& q2, double& u,
        double &ub, double &d, double &db, double &s, double &sb, double &c,
        double &b, double &g);*/
/*    void fdss_(int& hadron, int& charge, int& order, double& z,
              double& scalesqr, int& parton, double& result);

    extern struct{
        double fini;
    } fragini_;
}*/

int main(int argc, char* argv[]){

    using namespace params;

    double zmin = str_to_double(argv[1]);
    double zmax = str_to_double(argv[2]);
    double zstep = str_to_double(argv[3]);
    col = argv[4];
    b = str_to_double(argv[5]);
    incoming = argv[6];
    outgoing = argv[7];
    string rc = argv[8];
    p = str_to_double(argv[9]);
    double muratio = str_to_double(argv[10]);
    mu2 = Sq(muratio*p);
    if(rc.compare("fixed") == 0) alpha_s_running = FIXED;
    else if(rc.compare("mom") == 0) alpha_s_running = MOM;
    else if(rc.compare("parent") == 0) alpha_s_running = PARENT;
    else if(rc.compare("daughter") == 0) alpha_s_running = DAUGHTER;
    else if(rc.compare("mixed") == 0) alpha_s_running = MIXED;
    else if(rc.compare("mixedbd") == 0) alpha_s_running = MIXEDBD;
    else if(rc.compare("smallest") == 0) alpha_s_running = SMALLEST;
    else cout << "Invalid running coupling choice." << endl;

    ifstream datafile(TAfile.c_str());
    if(!datafile.is_open()){
        cerr << "Error opening the TA file" << endl;
        exit(1);
    }
    while(!datafile.eof()){
        string line;
        getline(datafile, line);
        if(line.size() == 0){
            break;
        }
        int loc = 0;
        while(line.compare(loc,1,",") != 0){
            loc++;
        }
        double this_b = str_to_double(line.substr(0, loc));
        if(this_b == b){
            TA = str_to_double(line.substr(loc+1, line.size()-loc-1));
            break;
        }
    }
    datafile.close();

    if(incoming.compare("g") == 0){
        if(outgoing.compare("g") == 0) with_gg = true;
        else with_gq = true;
    }
    else{
        if(outgoing.compare("g") == 0) with_gl = true;
        else{
            with_Nc = true;
            with_CF = true;
        }
    }

    init();

    double z = zmax;
    double k = p/zmax;
    //double k = 20.1;
    while(z >= zmin - 0.00001){
    //while(k <= 50.00001){

      //mu2 = Sq(4*k);

      k = p/z;
      double xp=(k/SQRTS)*exp(yh);
      double xg=(k/SQRTS)*exp(-yh);
      init_interp(xp,xg,k);
      cout << z << "," << k << "," << sigma_LO_k(k,xp) << "," << sigma_NLO_k(k,xp) << endl;
      //cout << k << "," << sigma_LO_k(k,xp) << "," << sigma_NLO_k(k,xp) << endl;
      clear_interp();
      z-=zstep;
      //k += 0.1;
   }

    return 0;
}
