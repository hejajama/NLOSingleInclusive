#include "params.hpp"
#include "utils.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

namespace params{

double lookup_TA(double b){
  ifstream datafile(TAfile.c_str());
  if(!datafile.is_open()){
    cerr << "Error opening the TA file" << endl;
    exit(1);
  }
  double TA = 0;
  while(!datafile.eof()){
    string line;
    getline(datafile, line);
    if(line.size() == 0){
      break;
    }
    int loc = 0;
    while(line.compare(loc,1,",") != 0){
      loc++;
    }
    double this_b = str_to_double(line.substr(0, loc));
    if(this_b == b){
      TA = str_to_double(line.substr(loc+1, line.size()-loc-1));
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

  return rp;
}

}
