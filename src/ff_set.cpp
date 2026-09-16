#include "ff_set.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;
using namespace params;

namespace {
  // Splits on literal commas, no whitespace trimming -- set names never
  // contain spaces, so "a,b" -> {"a","b"} and "a" -> {"a"}.
  vector<string> split_comma(const string& s){
    vector<string> parts;
    size_t start = 0;
    while(true){
      size_t comma = s.find(',', start);
      if(comma == string::npos){
        parts.push_back(s.substr(start));
        break;
      }
      parts.push_back(s.substr(start, comma-start));
      start = comma+1;
    }
    return parts;
  }
}


FfSet::FfSet(const string& name): name_(name){
  LHAPDF::setVerbosity(0);
  vector<string> parts = split_comma(name);
  if(parts.size() == 1){
    ff_.reset(LHAPDF::mkPDF(parts[0], 0));
  }
  else if(parts.size() == 2){
    // Two-set averaging mode (see ff_set.hpp) -- e.g. a pi0 FF approximated
    // as (D_pi+ + D_pi-)/2.
    ff_.reset(LHAPDF::mkPDF(parts[0], 0));
    ff2_.reset(LHAPDF::mkPDF(parts[1], 0));
  }
  else{
    cerr << "Bad FF set name (expected one name, or two separated by a "
            "single comma): " << name << endl;
    exit(1);
  }
}


double FfSet::eval(const LHAPDF::PDF& ff, const RunParameters& rp, double z, double Q2) const{
    // 1: d, -1: dbar
    // 2: u, -2: ubar
    // 3: s, -3: sbar
    // (mirrors PdfSet::xf's flavor mapping; gluon is 21 here rather than
    // PdfSet::xf's 0 -- both are LHAPDF aliases for the gluon, this just
    // keeps the value this branch always returned)
    const string &outgoing = rp.outgoing;
    if(outgoing.compare("g") == 0) return ff.xfxQ2(21,z,Q2);
    if(outgoing.compare("d") == 0) return ff.xfxQ2(1,z,Q2);
    if(outgoing.compare("u") == 0) return ff.xfxQ2(2,z,Q2);
    if(outgoing.compare("s") == 0) return ff.xfxQ2(3,z,Q2);
    if(outgoing.compare("dbar") == 0) return ff.xfxQ2(-1,z,Q2);
    if(outgoing.compare("ubar") == 0) return ff.xfxQ2(-2,z,Q2);
    if(outgoing.compare("sbar") == 0) return ff.xfxQ2(-3,z,Q2);
    if(outgoing.compare("q") == 0){
        // Sum over all Nf light quark AND antiquark flavors (params::Nf),
        // mirroring PdfSet::xf's rp.incoming=="q" handling. Only valid to
        // combine with a *single* explicit incoming flavor (or
        // rp.incoming=="g"/rp.channel!=QQ) -- see the rp.incoming=="q" &&
        // rp.outgoing=="q" rejection in sigma_hadron.cpp for why summing
        // both sides at once is wrong for the QQ channel.
        double sum = 0;
        for(int flavor=1; flavor<=Nf; flavor++){
            sum += ff.xfxQ2(flavor,z,Q2) + ff.xfxQ2(-flavor,z,Q2);
        }
        return sum;
    }

    cerr << "Bad rp.outgoing= " << outgoing << endl;
    exit(1);
}


double FfSet::zD(const RunParameters& rp, double z, double Q2) const{
    if(z<0 || z>1){
        cerr << "Bad z= " << z << endl;
        exit(1);
    }

    double result = eval(*ff_, rp, z, Q2);
    if(ff2_) result = 0.5*(result + eval(*ff2_, rp, z, Q2));
    return result;
}
