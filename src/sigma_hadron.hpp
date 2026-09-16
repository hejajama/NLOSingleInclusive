#pragma once

#include "dipole_amplitude.hpp"
#include "ff_set.hpp"
#include "params.hpp"
#include "pdf_set.hpp"

// Hadron-level single-inclusive cross section (arXiv:2310.06640 sec. V):
// convolves the parton-level LO/NLO cross sections (sigma_LO.hpp/
// sigma_NLO.hpp) with a collinear fragmentation function D(z,Q^2)
// (ff_set.hpp) over the momentum fraction z = p_h/k, where p_h is the
// observed hadron's transverse momentum and k is the parton's:
//
//   sigma_LO/NLO(p_h) = \int_{zmin}^{1} dz/z^2 * ff.zD(rp,z,rp.mu2_ff)
//                          * sigma_LO/NLO_k(rp, pdf, ..., k=p_h/z, xp(z))
// 
// xp(z)/xg(z) are recomputed at each z from k=p_h/z the same way main.cpp's
// parton-level loop computes them from k=p/z (Eqs. 5-6), which is why
// sqrts/y (not otherwise part of RunParameters) are taken as explicit
// arguments here. zmin is the lower cutoff of the (otherwise continuous)
// z integration -- its upper bound is always 1.
//
// The z integral itself is a fixed n-point Gauss-Legendre rule (GSL's
// gsl_integration_glfixed), not an adaptive one: each node costs a full
// NLO PointTables build (point_tables.cpp's per-r xi-convolution and
// per-(r,y) NLO-coefficient tables -- the same cost a single parton-level
// NLO point has), so an adaptive rule's node count -- and hence runtime --
// isn't something the caller can bound in advance. n directly trades
// accuracy for a predictable runtime of n times that per-point cost; see
// cli.hpp's --z-points and the validation note in sigma_hadron.cpp for how
// its default was chosen. LO is cheap (DipoleAmplitude::S(r,Y) directly,
// no PointTables) and shares the same n z-nodes as NLO, one PointTables
// build serving both per node.
struct HadronSigma{
  double LO;
  double NLO;
};

// Returns dN / d^2 p_h dy. To get the cross section, in pp collisions  this should be multiplied by params::sigma0/2
// in pA, one should instead perform the d^2b integral (which in pp case is replaced by \int d^2b -> simga_0/2)
HadronSigma sigma_hadron_ph(const params::RunParameters& rp, const PdfSet& pdf, const FfSet& ff,
                             const DipoleAmplitude& dipole, double p_h, double zmin,
                             double sqrts, double y, int n_zpoints);
