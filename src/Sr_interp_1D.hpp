#pragma once

// Cubic spline in r of S(r,Y) at a single, fixed rapidity Y — a cheaper
// slice of Sr_interp_2D used inside the r-integrals that run at fixed Y
// (e.g. the NLO coefficient functions in func.cpp).
double Sr_interp_1D(double r);

void init_Sr_interp_1D(double y);
void clear_Sr_interp_1D();
