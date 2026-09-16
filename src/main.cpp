#include "cli.hpp"
#include "common.hpp"
#include "dipole_amplitude.hpp"
#include "ff_set.hpp"
#include "params.hpp"
#include "pdf_set.hpp"
#include "point_tables.hpp"
#include "sigma_LO.hpp"
#include "sigma_NLO.hpp"
#include "sigma_hadron.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace params;

// CLI: see cli.hpp. Named flags (--zmin --col --incoming --outgoing --rc
// --pt --muratio, plus --zmax/--zstep unless --level hadron, plus --b if
// --col pA, plus optional --sqrts/--y/--bk-proton/--bk-nucleus/--sigma02)
// are preferred; the original positional form (zmin zmax zstep col b
// incoming outgoing rc p muratio) is still accepted for backward
// compatibility.
//
// Scans the parton-level LO+NLO single-inclusive cross section
// (arXiv:2310.06640 sec. IV) over z = p/k (p is the fixed CLI argument;
// k is the parton transverse momentum, p = z*k per the paragraph after
// Eq. 7). This is --level parton (the default, see cli.hpp) -- z is swept
// directly rather than integrated over, so it always reports parton-level
// results at k = p/z.
//
// --level hadron instead convolves that same parton-level cross section
// with a fragmentation function (sec. V) over z, reporting a single
// hadron-level result at the fixed hadron transverse momentum p_h = --pt
// (see sigma_hadron.hpp).
int main(int argc, char* argv[]){

    const cli::Args args = cli::parse(argc, argv);

    double mu2 = Sq(args.muratio*args.pt);
    running_types alpha_s_running = parse_alpha_s_running(args.rc);

    const RunParameters rp = make_run_parameters(args.col, args.b, args.pt, args.incoming, args.outgoing,
                                                  alpha_s_running, mu2,
                                                  args.bk_proton, args.bk_nucleus, args.sigma0);

    PdfSet pdf(args.pdf_set);
    DipoleAmplitude dipole(rp);
    gsl_set_error_handler(&gsl_error_handler);

    if(args.level == "parton"){
      double z = args.zmax;
      double k = args.pt/args.zmax;
     cout << "z,k [GeV],dN_LO/d^2kdy [1/GeV^2],dN_NLO/d^2kdy [1/GeV^2]" << endl;
      while(z >= args.zmin - 0.00001){

        k = args.pt/z;
        double xp=(k/args.sqrts)*exp(args.y);   // Eq. 6
        double xg=(k/args.sqrts)*exp(-args.y);  // Eq. 5 (Xg in the paper's notation)
        PointTables tables(rp, pdf, dipole, xp, xg, k);
        cout << z << "," << k << "," << sigma_LO_k(rp,pdf,dipole,xg,k,xp) << "," << sigma_NLO_k(rp,pdf,tables,k,xp) << endl;
        z-=args.zstep;
      }
    }
    else{  // "hadron"
        FfSet ff(args.ff_set);
        double p_h = args.pt;
        HadronSigma sigma = sigma_hadron_ph(rp, pdf, ff, dipole, p_h, args.zmin, args.sqrts, args.y, args.z_points);
        // in pp case, multiply by sigma0/2 to get the cross section 
        if (args.col == "pp") {
            sigma.LO *= rp.sigma0/2.0;
            sigma.NLO *= rp.sigma0/2.0;
            cout << "pT [GeV],d sigma_LO / d^2p_Tdy [1/GeV^4],d sigma_NLO / d^2p_Tdy [1/GeV^4]" << endl;
        }
        else 
        {
            cout << "pT [GeV],dN_LO / d^2p_Tdy [1/GeV^2],dN_NLO / d^2p_Tdy [1/GeV^2]" << endl;
        }

      cout << p_h << "," << sigma.LO << "," << sigma.NLO << endl;
    }

    return 0;
}
