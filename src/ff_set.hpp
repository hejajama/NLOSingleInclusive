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
  // name is either a single LHAPDF FF set name, or two names separated by
  // a comma (no spaces), e.g. "NNFF10_PIp_nlo,NNFF10_Pim_nlo" -- zD() then
  // returns the average of the two sets' z*D(z,Q^2). LHAPDF has no
  // charge-neutral pion (pi0) fragmentation function, so this is how a
  // pi0 FF is approximated: isospin symmetry gives
  // D_pi0 = (D_pi+ + D_pi-)/2. See README.md, "pi0 fragmentation
  // functions".
  explicit FfSet(const std::string& name);
  FfSet(const FfSet&) = delete;
  FfSet& operator=(const FfSet&) = delete;

  // z*D(z,Q^2) for rp.outgoing fragmenting into the set's hadron.
  // rp.outgoing=="g" reads the gluon fragmentation function directly;
  // rp.outgoing can also name a specific light quark/antiquark flavor
  // ("u"/"d"/"s"/"ubar"/"dbar"/"sbar", mirroring PdfSet::xf's rp.incoming
  // flavor tokens) to read just that flavor's FF; rp.outgoing=="q" sums
  // the fragmentation functions of all Nf light quark AND antiquark
  // flavors (params::Nf) instead. For the QQ channel (see Channel in
  // params.hpp), rp.incoming=="q" && rp.outgoing=="q" together is wrong
  // (it mixes flavors that shouldn't mix -- see sigma_hadron.cpp) and is
  // rejected by sigma_hadron_ph; use matching explicit flavor pairs
  // (--incoming u --outgoing u, etc.) instead and sum the results -- see
  // README.md, "Hadron-level flavor sums".
  double zD(const params::RunParameters& rp, double z, double Q2) const;

  // The name (or comma-separated pair of names) this instance was built
  // from -- lets a caller build an independent FfSet wrapping the same
  // set(s) (e.g. one per thread, since LHAPDF::PDF is not safe to share
  // across threads: see sigma_hadron.cpp). Round-trips through the FfSet
  // constructor unchanged, averaging mode included.
  const std::string& name() const { return name_; }

private:
  double eval(const LHAPDF::PDF& ff, const params::RunParameters& rp, double z, double Q2) const;

  std::string name_;
  std::unique_ptr<LHAPDF::PDF> ff_;
  std::unique_ptr<LHAPDF::PDF> ff2_;  // non-null only in two-set averaging mode
};
