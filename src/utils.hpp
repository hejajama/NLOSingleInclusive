#pragma once

#include <string>
#include <LHAPDF/LHAPDF.h>

// Collinear PDF set used for the incoming/outgoing parton (params::incoming).
namespace pdf{
  extern LHAPDF::PDF* lhpdf;
}

// PDF x*f(x,Q^2) for params::incoming, from LHAPDF.
double xf(double x, double Q2);

// Position-space running coupling alpha_s(r) (Henri's parametrization).
double alpha_s_pos(double r);

// Momentum-space running coupling alpha_s(Q^2), one-loop.
double alpha_s_mom(double Q2);

double str_to_double(std::string str);
int str_to_int(std::string str);

// GSL error handler installed in init(): ignores a few benign GSL warnings
// (max-iterations/tolerance/underflow/rounding) and aborts on anything else.
void gsl_error_handler(const char * reason,
                       const char * file,
                       int line,
                       int error);
