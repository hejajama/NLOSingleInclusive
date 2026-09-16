#include "pdf_set.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;
using namespace params;

namespace {
  // Splits on literal '+', no whitespace trimming -- flavor tokens never
  // contain spaces, so "u+ubar" -> {"u","ubar"} and "u" -> {"u"}.
  vector<string> split_plus(const string& s){
    vector<string> parts;
    size_t start = 0;
    while(true){
      size_t plus = s.find('+', start);
      if(plus == string::npos){
        parts.push_back(s.substr(start));
        break;
      }
      parts.push_back(s.substr(start, plus-start));
      start = plus+1;
    }
    return parts;
  }

  // x*f(x,Q^2) for a single named flavor token (never "q" -- that's
  // PdfSet::xf's own all-Nf-flavor aggregate, handled separately).
  double single_flavor_xf(const LHAPDF::PDF& pdf, const string& token, double x, double Q2){
    // 0: gluon
    // 1: d, -1: dbar
    // 2: u, -2: ubar
    // 3: s, -3: sbar
    if(token.compare("g") == 0) return pdf.xfxQ2(0,x,Q2);
    if(token.compare("d") == 0) return pdf.xfxQ2(1,x,Q2);
    if(token.compare("u") == 0) return pdf.xfxQ2(2,x,Q2);
    if(token.compare("s") == 0) return pdf.xfxQ2(3,x,Q2);
    if(token.compare("dbar") == 0) return pdf.xfxQ2(-1,x,Q2);
    if(token.compare("ubar") == 0) return pdf.xfxQ2(-2,x,Q2);
    if(token.compare("sbar") == 0) return pdf.xfxQ2(-3,x,Q2);

    cerr << "Bad rp.incoming flavor token: " << token << endl;
    exit(1);
  }
}


PdfSet::PdfSet(const string& name): name_(name){
  LHAPDF::setVerbosity(0);
  pdf_.reset(LHAPDF::mkPDF(name, 0));
}


double PdfSet::xf(const RunParameters& rp, double x, double Q2) const{
    if(x<0 || x>1){
        cerr << "Bad x= " << x << endl;
        exit(1);
    }

    const string &incoming = rp.incoming;
    if(incoming.compare("q") == 0){
        // Sum over all Nf light quark AND antiquark flavors (params::Nf),
        // mirroring FfSet::zD's rp.outgoing=="q" handling -- see
        // docs/pdf_evaluation_bug.md for why this branch didn't always
        // exist (it used to fall through with particle_flag uninitialized).
        double sum = 0;
        for(int flavor=1; flavor<=Nf; flavor++){
            sum += pdf_->xfxQ2(flavor,x,Q2) + pdf_->xfxQ2(-flavor,x,Q2);
        }
        return sum;
    }

    // A single flavor token ("u"), or several joined by '+' ("u+ubar") to
    // sum just those -- e.g. for a self-conjugate hadron (pi0, ...) where
    // the outgoing FF is the same for a flavor and its antiflavor
    // (D_u==D_ubar), summing the matching incoming pair here lets one
    // sigma_hadron_ph call produce both flavors' contributions while
    // building PointTables (the expensive part -- NLOCoefficients'
    // channel-only-dependent double integrals) only once instead of twice.
    // See README.md, "pi0 fragmentation functions", and the rejection of
    // multi-flavor incoming + rp.outgoing=="q" in sigma_hadron.cpp (that
    // combination is *not* safe in general -- only exploitable when the FF
    // really is flavor-independent across the summed set, which the
    // caller, not this function, is responsible for knowing).
    double sum = 0;
    for(const string& token : split_plus(incoming)){
      sum += single_flavor_xf(*pdf_, token, x, Q2);
    }
    return sum;
}
