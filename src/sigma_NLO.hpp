#pragma once

#include "params.hpp"
#include "point_tables.hpp"
#include "utils.hpp"

// NLO single-inclusive cross section correction (arXiv:2310.06640), obtained
// by Hankel-transforming sigma_NLO_r(r,k,xp), which combines the LO piece
// (via Sr_0) with the xi-convolved NLO correction (tables.xi_convolution).
double sigma_NLO_k(const params::RunParameters& rp, const PdfSet& pdf,
                    const PointTables& tables, double k, double xp);
