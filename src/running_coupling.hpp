#pragma once

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
