#pragma once

#include "params.hpp"
#include "pdf_set.hpp"
#include "point_tables.hpp"

// NLO single-inclusive PARTON-level cross section correction
// (arXiv:2310.06640 Eqs. 9 and 11a-11c, sec. IV -- same parton/hadron-level
// caveat as sigma_LO.hpp), obtained by Hankel-transforming
// sigma_NLO_r(r,k,xp), which combines the LO piece (via Sr_0, Eq. 13/14)
// with the xi-convolved NLO correction (tables.xi_convolution).
double sigma_NLO_k(const params::RunParameters& rp, const PdfSet& pdf,
                    const PointTables& tables, double k, double xp);
