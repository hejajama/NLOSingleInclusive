#pragma once

#include "params.hpp"
#include "pdf_set.hpp"
#include "point_tables.hpp"

// LO single-inclusive PARTON-level cross section (arXiv:2310.06640
// Eqs. 7a/7b, sec. IV -- the fragmentation-function convolution to the
// hadron level, sec. V, is not implemented here), as a function of
// transverse momentum k and parton momentum fraction xp. Obtained by
// Hankel-transforming the LO dipole-amplitude convolution sigma_LO_r(r,xp)
// (Eq. 7a for the quark channel, rp.with_Nc/with_CF; Eq. 7b for the gluon
// channel, rp.with_gg).
double sigma_LO_k(const params::RunParameters& rp, const PdfSet& pdf,
                   const PointTables& tables, double k, double xp);
