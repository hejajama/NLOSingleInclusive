#include "Sr_interp_1D.hpp"

#include <algorithm>

using namespace std;

DipoleAmplitude1DSlice::DipoleAmplitude1DSlice(const DipoleAmplitude& dipole, double Y)
  : minr_(dipole.min_r()), maxr_(dipole.max_r())
{
  int rpoints = dipole.r_points();
  const vector<double>& rvals = dipole.r_values();
  vector<double> Srvals(rpoints);
  for(int i=0; i<rpoints; i++){
    double r = rvals[i];
    Srvals[i] = max(0., dipole.S(r, Y));
  }
  spline_.build(rvals.data(), Srvals.data(), rpoints);
}


double DipoleAmplitude1DSlice::operator()(double r) const{
  if(r<minr_) return 1;
  if(r>maxr_) return 0;
  return spline_.eval(r);
}
