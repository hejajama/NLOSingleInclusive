#pragma once

#include "params.hpp"
#include "pdf_set.hpp"
#include "point_tables.hpp"

// NLO single-inclusive PARTON-level cross section correction
// (arXiv:2310.06640 Eqs. 9 and 11a-11c, sec. IV -- same parton/hadron-level
// caveat as sigma_LO.hpp), obtained by Hankel-transforming
// sigma_NLO_r(r,k,xp), which combines the LO piece (via Sr_0, Eq. 13/14,
// evaluated at the paper's X0 -- see sigma_LO_k_X0's doc comment) with the
// xi-convolved NLO correction (tables.xi_convolution). include_lo_baseline
// defaults to true, keeping this the full total (Eq. 7+9/11) for existing
// callers; pass false to get only the xi-convolution correction term
// (Eq. 9/11 alone), e.g. when the caller (sigma_hadron_ph) wants to add
// the X0 baseline itself, via sigma_LO_k_X0, instead of getting it bundled
// in here.
double sigma_NLO_k(const params::RunParameters& rp, const PdfSet& pdf,
                    const PointTables& tables, double k, double xp,
                    bool include_lo_baseline = true);
