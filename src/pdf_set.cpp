#include "pdf_set.hpp"

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace params;

PdfSet::PdfSet(const string& name){
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
    int particle_flag;
    const string &incoming = rp.incoming;
    if(incoming.compare("g") == 0) particle_flag = 0;
    else if(incoming.compare("d") == 0) particle_flag = 1;
    else if(incoming.compare("u") == 0) particle_flag = 2;
    else if(incoming.compare("s") == 0) particle_flag = 3;
    else if(incoming.compare("dbar") == 0) particle_flag = -1;
    else if(incoming.compare("ubar") == 0) particle_flag = -2;
    else if(incoming.compare("sbar") == 0) particle_flag = -3;

    return pdf_->xfxQ2(particle_flag,x,Q2);
}
