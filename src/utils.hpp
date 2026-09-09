#pragma once

#include "params.hpp"

#include <string>
#include <LHAPDF/LHAPDF.h>

// Collinear PDF set used for the incoming/outgoing parton (rp.incoming).
namespace pdf{
  extern LHAPDF::PDF* lhpdf;
}

// PDF x*f(x,Q^2) for rp.incoming, from LHAPDF.
double xf(const params::RunParameters& rp, double x, double Q2);

// Position-space running coupling alpha_s(r) (Henri's parametrization).
// Depends only on fixed physics constants (params::beta0/LambdaQCD), not on
// the run configuration.
double alpha_s_pos(double r);

// Momentum-space running coupling alpha_s(Q^2), one-loop. Same as above:
// no run-configuration dependence.
double alpha_s_mom(double Q2);

double str_to_double(std::string str);
int str_to_int(std::string str);

// GSL error handler installed in init(): ignores a few benign GSL warnings
// (max-iterations/tolerance/underflow/rounding) and aborts on anything else.
void gsl_error_handler(const char * reason,
                       const char * file,
                       int line,
                       int error);
