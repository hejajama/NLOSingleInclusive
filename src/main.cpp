#include "cli.hpp"
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

// CLI: see cli.hpp. Named flags (--zmin --zmax --zstep --col --b --incoming
// --outgoing --rc --p --muratio, plus optional --sqrts/--y/--bk-proton/
// --bk-nucleus/--sigma02) are preferred; the original positional form
// (zmin zmax zstep col b incoming outgoing rc p muratio) is still
// accepted for backward compatibility.
//
// Scans the parton-level LO+NLO single-inclusive cross section
// (arXiv:2310.06640 sec. IV) over z = p/k (p is the fixed CLI argument;
// k is the parton transverse momentum, p = z*k per the paragraph after
// Eq. 7). No fragmentation-function convolution to the hadron level
// (sec. V) is performed -- z is swept directly rather than integrated
// over, so this always reports parton-level results at k = p/z.
int main(int argc, char* argv[]){

    const cli::Args args = cli::parse(argc, argv);

    double mu2 = Sq(args.muratio*args.p);
    running_types alpha_s_running = parse_alpha_s_running(args.rc);

    const RunParameters rp = make_run_parameters(args.col, args.b, args.p, args.incoming, args.outgoing,
                                                  alpha_s_running, mu2,
                                                  args.bk_proton, args.bk_nucleus, args.sigma0);

    PdfSet pdf(pdfname);
    DipoleAmplitude dipole(rp);
    gsl_set_error_handler(&gsl_error_handler);

    double z = args.zmax;
    double k = args.p/args.zmax;
    while(z >= args.zmin - 0.00001){

      k = args.p/z;
      double xp=(k/args.sqrts)*exp(args.y);   // Eq. 6
      double xg=(k/args.sqrts)*exp(-args.y);  // Eq. 5 (Xg in the paper's notation)
      PointTables tables(rp, pdf, dipole, xp, xg, k);
      cout << z << "," << k << "," << sigma_LO_k(rp,pdf,tables,k,xp) << "," << sigma_NLO_k(rp,pdf,tables,k,xp) << endl;
      z-=args.zstep;
   }

    return 0;
}
