#include "pdf_set.hpp"

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace params;

PdfSet::PdfSet(const string& name): name_(name){
  LHAPDF::setVerbosity(0);
  pdf_.reset(LHAPDF::mkPDF(name, 0));
}


double PdfSet::xf(const RunParameters& rp, double x, double Q2) const{
    if(x<0 || x>1){
        cerr << "Bad x= " << x << endl;
        exit(1);
    }

    // 0: gluon
    // 1: d, -1: dbar
    // 2: u, -2: ubar
    // 3: s, -3: sbar
    const string &incoming = rp.incoming;
    if(incoming.compare("g") == 0) return pdf_->xfxQ2(0,x,Q2);
    if(incoming.compare("d") == 0) return pdf_->xfxQ2(1,x,Q2);
    if(incoming.compare("u") == 0) return pdf_->xfxQ2(2,x,Q2);
    if(incoming.compare("s") == 0) return pdf_->xfxQ2(3,x,Q2);
    if(incoming.compare("dbar") == 0) return pdf_->xfxQ2(-1,x,Q2);
    if(incoming.compare("ubar") == 0) return pdf_->xfxQ2(-2,x,Q2);
    if(incoming.compare("sbar") == 0) return pdf_->xfxQ2(-3,x,Q2);
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

    cerr << "Bad rp.incoming= " << incoming << endl;
    exit(1);
}
