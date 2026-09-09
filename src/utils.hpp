#pragma once

#include "params.hpp"

#include <memory>
#include <string>
#include <LHAPDF/LHAPDF.h>

// The collinear PDF set used for the incoming/outgoing parton
// (rp.incoming). Built once per run in main() from params::pdfname.
class PdfSet{
public:
  explicit PdfSet(const std::string& name);
  PdfSet(const PdfSet&) = delete;
  PdfSet& operator=(const PdfSet&) = delete;

  // x*f(x,Q^2) for rp.incoming.
  double xf(const params::RunParameters& rp, double x, double Q2) const;

private:
  std::unique_ptr<LHAPDF::PDF> pdf_;
};

// Position-space running coupling alpha_s(r) (Henri's parametrization).
// Depends only on fixed physics constants (params::beta0/LambdaQCD), not on
// the run configuration.
double alpha_s_pos(double r);

// Momentum-space running coupling alpha_s(Q^2), one-loop. Same as above:
// no run-configuration dependence.
double alpha_s_mom(double Q2);

double str_to_double(std::string str);
int str_to_int(std::string str);

// GSL error handler installed in main(): ignores a few benign GSL warnings
// (max-iterations/tolerance/underflow/rounding) and aborts on anything else.
void gsl_error_handler(const char * reason,
                       const char * file,
                       int line,
                       int error);
