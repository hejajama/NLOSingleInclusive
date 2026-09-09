#pragma once

#include <vector>

// The BK-evolved dipole amplitude S(r,Y), read from a solved-BK-equation
// grid file (arXiv:2310.06640 sec. 2). "amplitude" holds that grid plus the
// r/Y sampling it was tabulated on; init_bksol() populates it by reading
// params::bksolpp / params::bksolpA.
namespace amplitude{
  extern double minr, maxr, r_mult, x0, Y0threshold;
  extern int rpoints, ypoints;
  extern double *yvals, *rvals;
  extern std::vector<std::vector<double> > Sr;  // Sr[Y][r]=S(Y,r)
}

// Analytic (GBW-like) initial-condition dipole amplitude S(r) at the BK
// evolution starting rapidity, used directly by the NLO cross section
// (sigma_NLO.cpp) rather than through the interpolated grid.
double Sr_0(double r);

// Loads amplitude::Sr (and its r/Y grid) from the BK solution file selected
// by params::col/params::b, and derives params::minlnr/maxlnr from it.
void init_bksol();
