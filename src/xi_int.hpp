#pragma once

#include "params.hpp"

// The xi (momentum-fraction) convolution of the NLO coefficient functions
// with the collinear PDF, arXiv:2310.06640 sec. 3 — this is what turns the
// per-channel coefficient functions into the NLO correction to sigma(r).
// xi_int_interp(r) is a cubic spline in r of that convolution integral at
// the current (xp, xg, k), built by init_xi_int_interp().
double xi_int_interp(double r);

void init_xi_int_interp(const params::RunParameters& rp, double xp, double xg, double k);
void clear_xi_int_interp();
