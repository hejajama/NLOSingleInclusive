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

  // x*f(x,Q^2) for rp.incoming. rp.incoming=="q" sums the PDFs of all Nf
  // light quark AND antiquark flavors (params::Nf) -- see
  // docs/pdf_evaluation_bug.md for the history of this branch.
  double xf(const params::RunParameters& rp, double x, double Q2) const;

  // The LHAPDF set name this instance was built from -- lets a caller build
  // an independent PdfSet wrapping the same set (e.g. one per thread, since
  // LHAPDF::PDF is not safe to share across threads: see sigma_hadron.cpp).
  const std::string& name() const { return name_; }

private:
  std::string name_;
  std::unique_ptr<LHAPDF::PDF> pdf_;
};
