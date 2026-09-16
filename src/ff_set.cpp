#include "ff_set.hpp"

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace params;

FfSet::FfSet(const string& name): name_(name){
  LHAPDF::setVerbosity(0);
  ff_.reset(LHAPDF::mkPDF(name, 0));
}


double FfSet::zD(const RunParameters& rp, double z, double Q2) const{
    if(z<0 || z>1){
        cerr << "Bad z= " << z << endl;
        exit(1);
    }

    const string &outgoing = rp.outgoing;
    if(outgoing.compare("g") == 0){
        return ff_->xfxQ2(21,z,Q2);
    }
    else if(outgoing.compare("q") == 0){
        double sum = 0;
        for(int flavor=1; flavor<=Nf; flavor++){
            sum += ff_->xfxQ2(flavor,z,Q2) + ff_->xfxQ2(-flavor,z,Q2);
        }
        return sum;
    }

    cerr << "Bad rp.outgoing= " << outgoing << endl;
    exit(1);
}
