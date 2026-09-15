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
//                          * sigma_LO/NLO_k(rp, pdf, tables(z), k=p_h/z, xp(z))
//
// xp(z)/xg(z) are recomputed at each z from k=p_h/z the same way main.cpp's
// parton-level loop computes them from k=p/z (Eqs. 5-6), which is why
// sqrts/y (not otherwise part of RunParameters) are taken as explicit
// arguments here. zmin is the lower cutoff of the (otherwise continuous)
// z integration -- its upper bound is always 1.
//
// LO and NLO are integrated separately, but LO's z-integral never builds a
// PointTables at all -- it only ever needs DipoleAmplitude::S(r,Y)
// (sigma_LO.cpp), never PointTables' NLO coefficient/xi-convolution tables
// (point_tables.cpp's per-r xi-convolution and per-(r,y) NLO-coefficient
// tables are what make building one expensive). See sigma_hadron.cpp.
struct HadronSigma{
  double LO;
  double NLO;
};

HadronSigma sigma_hadron_ph(const params::RunParameters& rp, const PdfSet& pdf, const FfSet& ff,
                             const DipoleAmplitude& dipole, double p_h, double zmin,
                             double sqrts, double y);
