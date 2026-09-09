#pragma once

#include "params.hpp"

// Bicubic (r,Y) interpolations of the func.hpp coefficient functions,
// built once per event by init_xi_interp() so the xi-convolution
// (xi_int.cpp) can evaluate them cheaply inside its xi integral.
double I2_interp(double r, double y);
double J_interp(double r, double y);
double J1_interp(double r, double y);
double H2_interp(double r, double y);
double H3_interp(double r, double y);
double H5_interp(double r, double y);
double K3_interp(double r, double y);
double Jv_interp(double r, double y);
double Jv2_interp(double r, double y);
double JJv_xi1_interp(double r, double y);

// xg: the gluon momentum fraction at the current kinematic point, used to
// bound the rapidity range the coefficient functions are tabulated over.
// Which coefficient functions actually get computed (vs. left at zero)
// depends on rp.with_CF/with_Nc/with_gl/with_gq/with_gg.
void init_xi_interp(const params::RunParameters& rp, double xg);
void clear_xi_interp();
