#pragma once

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
void init_xi_interp(double xg);
void clear_xi_interp();
