#pragma once

#include "params.hpp"

#include <gsl/gsl_integration.h>

// GSL integration workspaces used by func()'s phi/x double integral.
// Exposed here only because init_xi_interp() (xi_interp.cpp) allocates and
// frees them itself, once per rapidity slice, to amortize the cost across
// the r-loop at that slice — a cross-file lifetime dependency that a later
// refactor stage should own via a proper RAII type instead.
namespace func_tmp{
  extern gsl_integration_workspace *w_x, *w_phi;
#pragma omp threadprivate(w_x,w_phi)
}

// The nine NLO impact-factor coefficient functions of arXiv:2310.06640:
// I2, J1, H2, K3 (real corrections, various qq/qg/gq/gg channels), J, Jv,
// Jv2 (Nc-channel real/virtual pair), H3, H5 (gg-channel real corrections).
// Each numerically integrates func()'s underlying x,phi double integral
// (see func.cpp) with a channel-specific integrand selected by an internal
// flag. JJv_xi1 is the xi->1 subtraction term used only when
// params::with_xi1 is enabled. All depend on rp.alpha_s_running (the
// running-coupling prescription).
double I2(const params::RunParameters& rp, double r, double xi);
double J(const params::RunParameters& rp, double r, double xi);
double J1(const params::RunParameters& rp, double r, double xi);
double H2(const params::RunParameters& rp, double r, double xi);
double H3(const params::RunParameters& rp, double r, double xi);
double H5(const params::RunParameters& rp, double r, double xi);
double K3(const params::RunParameters& rp, double r, double xi);
double Jv(const params::RunParameters& rp, double r, double xi);
double Jv2(const params::RunParameters& rp, double r, double xi);
double JJv_xi1(const params::RunParameters& rp, double r);
