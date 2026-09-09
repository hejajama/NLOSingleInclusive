#pragma once

#include "params.hpp"
#include "pdf_set.hpp"
#include "point_tables.hpp"

// LO single-inclusive cross section (as a function of transverse momentum
// k and parton momentum fraction xp), arXiv:2310.06640, obtained by
// Hankel-transforming the LO dipole-amplitude convolution sigma_LO_r(r,xp).
double sigma_LO_k(const params::RunParameters& rp, const PdfSet& pdf,
                   const PointTables& tables, double k, double xp);
