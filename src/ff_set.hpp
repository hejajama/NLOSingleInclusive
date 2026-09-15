#pragma once

#include "params.hpp"

#include <memory>
#include <string>
#include <LHAPDF/LHAPDF.h>

// The collinear fragmentation-function set used to convolve the outgoing
// parton (rp.outgoing) into an observed hadron (sec. V, arXiv:2310.06640).
// LHAPDF ships fragmentation functions in the same grid format as PDFs, so
// this mirrors PdfSet (pdf_set.hpp) almost exactly -- built once per run in
// main() from a fragmentation-function set name such as "NNFF10_PIsum_nlo".
class FfSet{
public:
  explicit FfSet(const std::string& name);
  FfSet(const FfSet&) = delete;
  FfSet& operator=(const FfSet&) = delete;

  // z*D(z,Q^2) for rp.outgoing fragmenting into the set's hadron.
  // rp.outgoing=="g" reads the gluon fragmentation function directly;
  // rp.outgoing=="q" sums the fragmentation functions of all Nf light
  // quark AND antiquark flavors (params::Nf), since -- unlike rp.incoming
  // on the PdfSet side -- rp.outgoing never names a specific flavor (see
  // cli.hpp/params.cpp: outgoing is just "g" or "q").
  double zD(const params::RunParameters& rp, double z, double Q2) const;

private:
  std::unique_ptr<LHAPDF::PDF> ff_;
};
