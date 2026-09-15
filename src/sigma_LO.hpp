#pragma once

#include "dipole_amplitude.hpp"
#include "params.hpp"
#include "pdf_set.hpp"

// LO single-inclusive cross section (arXiv:2310.06640 Eqs. 7a/7b, sec. IV),
// as a function of transverse momentum k and parton momentum fraction xp,
// at rapidity Y = log(1/xg). Obtained by Hankel-transforming the LO
// dipole-amplitude convolution against DipoleAmplitude::S(r,Y) directly
// (Eq. 7a for the quark channel, rp.channel==Channel::QQ; Eq. 7b for the
// gluon channel, rp.channel==Channel::GG) -- no separate PointTables or
// DipoleAmplitude1DSlice needed: S(r,Y) is itself just a bicubic-spline
// lookup (dipole_amplitude.cpp), and LO never needs PointTables' NLO
// coefficient/xi-convolution tables.
double sigma_LO_k(const params::RunParameters& rp, const PdfSet& pdf,
                   const DipoleAmplitude& dipole, double xg, double k, double xp);
