#include "common.hpp"
#include "init.hpp"
#include "init_interp.hpp"
#include "params.hpp"
#include "sigma_LO.hpp"
#include "sigma_NLO.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
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

    double zmin = str_to_double(argv[1]);
    double zmax = str_to_double(argv[2]);
    double zstep = str_to_double(argv[3]);
    string col = argv[4];
    double b = str_to_double(argv[5]);
    string incoming = argv[6];
    string outgoing = argv[7];
    string rc = argv[8];
    double p = str_to_double(argv[9]);
    double muratio = str_to_double(argv[10]);
    double mu2 = Sq(muratio*p);
    running_types alpha_s_running = parse_alpha_s_running(rc);

    const RunParameters rp = make_run_parameters(col, b, p, incoming, outgoing,
                                                  alpha_s_running, mu2);

    init(rp);

    double z = zmax;
    double k = p/zmax;
    //double k = 20.1;
    while(z >= zmin - 0.00001){
    //while(k <= 50.00001){

      //mu2 = Sq(4*k);

      k = p/z;
      double xp=(k/SQRTS)*exp(yh);
      double xg=(k/SQRTS)*exp(-yh);
      init_interp(rp,xp,xg,k);
      cout << z << "," << k << "," << sigma_LO_k(rp,k,xp) << "," << sigma_NLO_k(rp,k,xp) << endl;
      //cout << k << "," << sigma_LO_k(rp,k,xp) << "," << sigma_NLO_k(rp,k,xp) << endl;
      clear_interp(rp);
      z-=zstep;
      //k += 0.1;
   }

    return 0;
}
