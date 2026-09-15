#include "dipole_amplitude.hpp"

#include "common.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <gsl/gsl_math.h>

using namespace std;
using namespace params;

namespace{
  // Finds `marker` in `line` and parses the first floating-point number
  // that follows it (skipping over "=", units labels, etc.) into `value`.
  // Used to pull the initial-condition parameters out of the BK solution
  // file's free-text header comment, e.g. marker "Q_s0^2" against
  // "Q_s0^2 = 0.0963203 GeV^2" -> value = 0.0963203. Returns false (leaving
  // `value` untouched) if the marker isn't in the line or no number follows
  // it.
  bool extract_after(const string& line, const string& marker, double& value){
    size_t pos = line.find(marker);
    if(pos == string::npos) return false;
    size_t start = pos + marker.size();
    while(start < line.size()
          && !std::isdigit(static_cast<unsigned char>(line[start]))
          && line[start] != '-' && line[start] != '.'){
      start++;
    }
    if(start >= line.size()) return false;
    try{
      value = std::stod(line.substr(start));
    } catch(const std::exception&){
      return false;
    }
    return true;
  }
}

InitialConditionParams read_initial_condition_header(const string& filename){
  InitialConditionParams icp;
  bool found_Qs02 = false, found_gamma = false, found_ec = false;

  ifstream datafile(filename.c_str());
  if(!datafile.is_open()){
    cerr << "Error opening the BK solution file " << filename
         << " to read its initial-condition header" << endl;
    exit(1);
  }
  int confid=0;
  string line;
  double val;
  while(confid<4 && getline(datafile,line)){
    if(extract_after(line, "Q_s0^2", val)){ icp.Qs02 = val; found_Qs02 = true; }
    if(extract_after(line, "\\gamma", val)){ icp.gamma = val; found_gamma = true; }
    if(extract_after(line, "coefficient of E inside Log is", val)){ icp.ec = val; found_ec = true; }
    if(extract_after(line, "Nf=", val)) icp.Nf = static_cast<int>(std::lround(val));
    if(line.substr(0,3)=="###") confid++;
  }
  datafile.close();

  // A file missing any of these in its header comment falls back to
  // params::default_Qs02/gamm/ec (the values Sr_0's initial condition used
  // to be hand-tuned to before this parser existed) rather than exiting --
  // but warns, since a silently wrong default is exactly the kind of
  // mismatch this parser was added to catch (see f75242d).
  if(!found_Qs02){
    cerr << "Warning: couldn't find Q_s0^2 in the initial-condition header of "
         << filename << " -- expected a line like \"# Initial condition: MV "
            "model, Q_s0^2 = ..., \\gamma = ..., coefficient of E inside Log "
            "is ..., x0=..., \\Lambda_QCD = ... GeV\"; using default Q_s0^2 = "
         << default_Qs02 << endl;
    icp.Qs02 = default_Qs02;
  }
  if(!found_gamma){
    cerr << "Warning: couldn't find \\gamma in the initial-condition header of "
         << filename << "; using default \\gamma = " << default_gamm << endl;
    icp.gamma = default_gamm;
  }
  if(!found_ec){
    cerr << "Warning: couldn't find \"coefficient of E inside Log\" in the "
            "initial-condition header of " << filename << "; using default "
            "value " << default_ec << endl;
    icp.ec = default_ec;
  }
  return icp;
}


double Sr_0(const RunParameters& rp, double r){
    if(rp.col.compare("pA") == 0){
        return exp(-0.125*rp.sigma0*Anucleus*rp.TA*pow(Sq(r)*rp.Qs02,rp.gamm)*log(1/(r*LambdaQCD)+rp.ec*M_E));
    }
    return exp(-0.25*pow(Sq(r)*rp.Qs02,rp.gamm)*log(1/(r*LambdaQCD)+rp.ec*M_E));
}


DipoleAmplitude::DipoleAmplitude(const RunParameters& rp){
  vector<vector<double>> Sr;
  vector<double> yvals_tmp;
  load_grid(rp, Sr, yvals_tmp);

  ypoints_ = yvals_tmp.size();
  yvals_ = yvals_tmp;

  rvals_.resize(rpoints_);
  for(int i=0; i<rpoints_; i++){
    rvals_[i] = minr_*gsl_pow_int(r_mult_, i);
  }
  maxr_ = rvals_[rpoints_-1];
  minlnr_ = log(0.5*minr_);
  maxlnr_ = log(2*maxr_);

  spline_.init(rpoints_, ypoints_);
  for(int i=0; i<ypoints_; i++){
    for(int j=0; j<rpoints_; j++){
      double tmp = max(0., Sr[i][j]);
      spline_.set(j, i, tmp);
    }
  }
  spline_.build(rvals_.data(), yvals_.data());
}


double DipoleAmplitude::S(double r, double Y) const{
  if(r<minr_) return 1;
  if(r>maxr_) return 0;
  // The bicubic spline can overshoot slightly (interpolation ripple) even
  // though S itself is a probability and must stay in [0,1] -- clamp rather
  // than let that leak into callers (e.g. sigma_LO.cpp squares this for the
  // gluon channel, where a tiny negative value is otherwise harmless but a
  // tiny >1 value is not).
  double s = spline_.eval(r, max(0.0, Y-std::log(1/x0_)));
  return clamp(s, 0.0, 1.0);
}


void DipoleAmplitude::load_grid(const RunParameters& rp,
                                vector<vector<double>>& Sr,
                                vector<double>& yvals_tmp){
  string bksol;
  if(rp.col.compare("pA") == 0){
    if(rp.b < 10){
      bksol = rp.bk_nucleus_prefix + to_string(rp.b).substr(0,1);
    }
    else{
      bksol = rp.bk_nucleus_prefix + to_string(rp.b).substr(0,2);
    }
  }
  else{
    bksol = rp.bk_proton_file;
  }

  InitialConditionParams icp = read_initial_condition_header(bksol);
  qs02_ = icp.Qs02;
  gamma_ = icp.gamma;
  ec_ = icp.ec;
  nf_ = icp.Nf;

  ifstream datafile(bksol.c_str());
  if(!datafile.is_open()){
    cerr << "Error opening the BK solution file" << endl;
    exit(1);
  }
  int confid=0;
  string line;
  while(confid<4 && getline(datafile,line)){
    if(line.substr(0,3)=="###"){
      switch(confid){
      case 0:
        minr_=std::stod(line.substr(3,line.length()-3));
        break;
      case 1:
        r_mult_=std::stod(line.substr(3,line.length()-3));
        break;
      case 2:
        rpoints_=std::stoi(line.substr(3,line.length()-3));
        break;
      case 3:
        x0_=std::stod(line.substr(3,line.length()-3));
        break;
      }
      confid++;
    }
  }
  vector<double> tmpvec;
  while(getline(datafile,line)){
    if(line.substr(0,3)=="###"){
      if(tmpvec.size()>0) Sr.push_back(tmpvec);
      yvals_tmp.push_back(std::stod(line.substr(3,line.length()-3)));
      tmpvec.clear();
      continue;
    }
    tmpvec.push_back(1-std::stod(line));
  }
  Sr.push_back(tmpvec);
  datafile.close();
}
