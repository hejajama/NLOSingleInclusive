#pragma once

#include "params.hpp"

#include <vector>

// The BK-evolved dipole amplitude S(r,Y), read from a solved-BK-equation
// grid file (arXiv:2310.06640 sec. 2). "amplitude" holds that grid plus the
// r/Y sampling it was tabulated on; init_bksol() populates it by reading
// params::RunParameters::col/b (via bksolpp/bksolpA).
//
// minlnr/maxlnr are the ln(r) integration bounds used by the NLO
// coefficient functions (func.cpp) -- derived from the grid's r-range, so
// they live here rather than in RunParameters.
namespace amplitude{
  extern double minr, maxr, r_mult, x0, Y0threshold;
  extern double minlnr, maxlnr;
  extern int rpoints, ypoints;
  extern double *yvals, *rvals;
  extern std::vector<std::vector<double> > Sr;  // Sr[Y][r]=S(Y,r)
}

// Analytic (GBW-like) initial-condition dipole amplitude S(r) at the BK
// evolution starting rapidity, used directly by the NLO cross section
// (sigma_NLO.cpp) rather than through the interpolated grid.
double Sr_0(const params::RunParameters& rp, double r);

// Loads amplitude::Sr (and its r/Y grid) from the BK solution file selected
// by rp.col/rp.b, and derives amplitude::minlnr/maxlnr from it.
void init_bksol(const params::RunParameters& rp);
