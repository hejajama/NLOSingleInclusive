#pragma once

// Bicubic interpolation of the BK-evolved dipole amplitude S(r,Y) over the
// full (r,Y) grid loaded by init_bksol() (see bksol_nlodisfit.hpp).
double Sr_interp_2D(double r, double y);

void init_Sr_interp_2D();
void clear_Sr_interp_2D();
