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
  if(rp.alpha_s_running==SMALLEST &&
     (rp.channel==Channel::QG || rp.channel==Channel::GQ || rp.channel==Channel::GG)){
    cerr << "Error: alpha_s_running=smallest is not implemented for the qg/gq/gg channels "
            "(the K1/K2/H2 coefficient functions never apply a running-coupling factor for "
            "this scheme). See docs/PAPER_MAPPING.md, \"Known coupling-scheme gaps\"." << endl;
    exit(1);
  }
  if(rp.alpha_s_running==DAUGHTER && rp.channel==Channel::QG){
    cerr << "Error: alpha_s_running=daughter is not implemented for the qg channel "
            "(the K1 coefficient function never applies a running-coupling factor for this "
            "scheme). See docs/PAPER_MAPPING.md, \"Known coupling-scheme gaps\"." << endl;
    exit(1);
  }
}


RunParameters make_run_parameters(string col, double b, double p,
                                   string incoming, string outgoing,
                                   running_types alpha_s_running, double mu2,
                                   const string& bk_proton, const string& bk_nucleus){
  RunParameters rp;
  rp.col = col;
  rp.b = b;
  rp.p = p;
  rp.incoming = incoming;
  rp.outgoing = outgoing;
  rp.alpha_s_running = alpha_s_running;
  rp.mu2 = mu2;
  rp.TA = lookup_TA(b);
  rp.bk_proton_file = bk_proton.empty() ? bksolpp : bk_proton;
  rp.bk_nucleus_prefix = bk_nucleus.empty() ? bksolpA : bk_nucleus;

  if(incoming.compare("g") == 0){
    rp.channel = (outgoing.compare("g") == 0) ? Channel::GG : Channel::GQ;
  }
  else{
    rp.channel = (outgoing.compare("g") == 0) ? Channel::QG : Channel::QQ;
  }

  validate_alpha_s_running(rp);

  return rp;
}

}
