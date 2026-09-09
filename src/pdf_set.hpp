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
