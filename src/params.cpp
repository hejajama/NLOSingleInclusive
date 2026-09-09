#include "params.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

namespace params{

double lookup_TA(double b){
  ifstream datafile(TAfile.c_str());
  if(!datafile.is_open()){
    cerr << "Error opening the TA file" << endl;
    exit(1);
  }
  double TA = 0;
  string line;
  while(getline(datafile, line)){
    if(line.size() == 0){
      break;
    }
    int loc = 0;
    while(line.compare(loc,1,",") != 0){
      loc++;
    }
    double this_b = std::stod(line.substr(0, loc));
    if(this_b == b){
      TA = std::stod(line.substr(loc+1, line.size()-loc-1));
      break;
    }
  }
  datafile.close();
  return TA;
}


running_types parse_alpha_s_running(const string& rc){
  if(rc.compare("fixed") == 0) return FIXED;
  else if(rc.compare("mom") == 0) return MOM;
  else if(rc.compare("parent") == 0) return PARENT;
  else if(rc.compare("daughter") == 0) return DAUGHTER;
  else if(rc.compare("mixed") == 0) return MIXED;
  else if(rc.compare("mixedbd") == 0) return MIXEDBD;
  else if(rc.compare("smallest") == 0) return SMALLEST;
  cout << "Invalid running coupling choice." << endl;
  return FIXED;
}


// Rejects (alpha_s_running, channel) combinations for which the
// corresponding NLO coefficient function (nlo_coefficients.cpp) never
// actually applies a running-coupling factor, in any of its own branches
// or the calling code's compensating logic (point_tables.cpp) -- see
// docs/PAPER_MAPPING.md, "Known coupling-scheme gaps", for the full
// derivation. Silently proceeding would compute those channels as if
// alpha_s were fixed to 1 for the affected term, which is physically
// wrong, not just imprecise -- so this exits rather than warns.
static void validate_alpha_s_running(const RunParameters& rp){
  if(rp.alpha_s_running==SMALLEST && (rp.with_gl || rp.with_gq || rp.with_gg)){
    cerr << "Error: alpha_s_running=smallest is not implemented for the qg/gq/gg channels "
            "(the J1/K3/H2 coefficient functions never apply a running-coupling factor for "
            "this scheme). See docs/PAPER_MAPPING.md, \"Known coupling-scheme gaps\"." << endl;
    exit(1);
  }
  if(rp.alpha_s_running==DAUGHTER && rp.with_gl){
    cerr << "Error: alpha_s_running=daughter is not implemented for the qg channel "
            "(the J1 coefficient function never applies a running-coupling factor for this "
            "scheme). See docs/PAPER_MAPPING.md, \"Known coupling-scheme gaps\"." << endl;
    exit(1);
  }
}


RunParameters make_run_parameters(string col, double b, double p,
                                   string incoming, string outgoing,
                                   running_types alpha_s_running, double mu2){
  RunParameters rp;
  rp.col = col;
  rp.b = b;
  rp.p = p;
  rp.incoming = incoming;
  rp.outgoing = outgoing;
  rp.alpha_s_running = alpha_s_running;
  rp.mu2 = mu2;
  rp.TA = lookup_TA(b);

  rp.with_Nc = false;
  rp.with_CF = false;
  rp.with_gl = false;
  rp.with_gq = false;
  rp.with_gg = false;
  if(incoming.compare("g") == 0){
    if(outgoing.compare("g") == 0) rp.with_gg = true;
    else rp.with_gq = true;
  }
  else{
    if(outgoing.compare("g") == 0) rp.with_gl = true;
    else{
      rp.with_Nc = true;
      rp.with_CF = true;
    }
  }

  validate_alpha_s_running(rp);

  return rp;
}

}
