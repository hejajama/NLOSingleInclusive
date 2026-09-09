#include "common.hpp"
#include "dipole_amplitude.hpp"
#include "params.hpp"
#include "pdf_set.hpp"
#include "point_tables.hpp"
#include "sigma_LO.hpp"
#include "sigma_NLO.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace params;

int main(int argc, char* argv[]){

    double zmin = std::stod(argv[1]);
    double zmax = std::stod(argv[2]);
    double zstep = std::stod(argv[3]);
    string col = argv[4];
    double b = std::stod(argv[5]);
    string incoming = argv[6];
    string outgoing = argv[7];
    string rc = argv[8];
    double p = std::stod(argv[9]);
    double muratio = std::stod(argv[10]);
    double mu2 = Sq(muratio*p);
    running_types alpha_s_running = parse_alpha_s_running(rc);

    const RunParameters rp = make_run_parameters(col, b, p, incoming, outgoing,
                                                  alpha_s_running, mu2);

    PdfSet pdf(pdfname);
    DipoleAmplitude dipole(rp);
    gsl_set_error_handler(&gsl_error_handler);

    double z = zmax;
    double k = p/zmax;
    while(z >= zmin - 0.00001){

      k = p/z;
      double xp=(k/SQRTS)*exp(yh);
      double xg=(k/SQRTS)*exp(-yh);
      PointTables tables(rp, pdf, dipole, xp, xg, k);
      cout << z << "," << k << "," << sigma_LO_k(rp,pdf,tables,k,xp) << "," << sigma_NLO_k(rp,pdf,tables,k,xp) << endl;
      z-=zstep;
   }

    return 0;
}
