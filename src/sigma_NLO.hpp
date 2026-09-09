#pragma once

// NLO single-inclusive cross section correction (arXiv:2310.06640), obtained
// by Hankel-transforming sigma_NLO_r(r,k,xp), which combines the LO piece
// (via Sr_0) with the xi-convolved NLO correction (xi_int_interp).
double sigma_NLO_k(double k, double xp);
