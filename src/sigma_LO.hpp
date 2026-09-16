#pragma once

#include "dipole_amplitude.hpp"
#include "params.hpp"
#include "pdf_set.hpp"

// LO single-inclusive cross section (arXiv:2310.06640 Eqs. 7a/7b, sec. IV),
// as a function of transverse momentum k and parton momentum fraction xp,
// at rapidity Y = log(1/xg) -- i.e. with the dipole BK-evolved to Xg,
// *not* the X0-frozen dipole the paper's own Eq. 7a/7b actually specifies
// (see sigma_LO_k_X0 below for that). Kept for reference/diagnostics.
// Obtained by Hankel-transforming the LO dipole-amplitude convolution
// against DipoleAmplitude::S(r,Y) directly (Eq. 7a for the quark channel,
// rp.channel==Channel::QQ; Eq. 7b for the gluon channel,
// rp.channel==Channel::GG) -- no separate PointTables or
// DipoleAmplitude1DSlice needed: S(r,Y) is itself just a bicubic-spline
// lookup (dipole_amplitude.cpp), and LO never needs PointTables' NLO
// coefficient/xi-convolution tables.
double sigma_LO_k(const params::RunParameters& rp, const PdfSet& pdf,
                   const DipoleAmplitude& dipole, double xg, double k, double xp);

// Eqs. 7a/7b evaluated at the dipole's initial-condition scale X0 (Eq.
// 13/14, via Sr_0()) rather than at Xg -- this is what the paper's own
// Eq. (7a)/(7b) actually calls for (the text right before those equations
// states dipole amplitudes in the LO term are evaluated "at the initial
// scale, X0, of the small-x evolution"; the unsubtracted-scheme rapidity
// evolution down to Xg only enters through the NLO xi-convolution's
// X(xi), Eq. 4). No DipoleAmplitude/BK grid needed at all, since Sr_0() is
// closed-form and X0-independent of k/xg.
double sigma_LO_k_X0(const params::RunParameters& rp, const PdfSet& pdf, double k, double xp);
