#pragma once

#include "params.hpp"

// Running-coupling prescriptions for alpha_s(r)/alpha_s(Q^2) used in the
// NLO impact factor (arXiv:2310.06640) -- which one applies to a given
// term is selected per params::RunParameters::alpha_s_running (see
// nlo_coefficients.cpp / point_tables.cpp / sigma_NLO.cpp). Depend only on
// fixed physics constants (params::beta0/LambdaQCD), not on the run
// configuration itself.

// Position-space running coupling alpha_s(r) (Henri's parametrization).
double alpha_s_pos(double r);

// Momentum-space running coupling alpha_s(Q^2), one-loop.
double alpha_s_mom(double Q2);

// The two alpha_s_running groupings that recur as an identical compound
// condition at several point_tables.cpp call sites (I1, H1, and the qq
// channel's C_F term): PARENT/DAUGHTER/SMALLEST all apply alpha_s_pos(r) to
// the same term there, and MIXED/MIXEDBD both apply a single momentum-space
// factor to it. These are NOT general "does scheme X apply a running-
// coupling factor here" predicates -- most call sites (every per-
// coefficient-function branch in nlo_coefficients.cpp, and the qg/gg
// blocks' lone PARENT checks in point_tables.cpp) give each scheme its own
// distinct treatment -- sometimes deliberately omitting one (see
// docs/PAPER_MAPPING.md, "Known coupling-scheme gaps") -- and are left as
// direct enum comparisons on purpose. Only use these two where a call
// site's condition is exactly this compound OR, unchanged; don't reach for
// them to "simplify" a single-value or differently-grouped check.
bool is_position_space_alpha_s(running_types rc);
bool is_mixed_alpha_s(running_types rc);
