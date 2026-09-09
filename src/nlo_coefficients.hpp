#pragma once

#include "Sr_interp_1D.hpp"
#include "dipole_amplitude.hpp"
#include "params.hpp"

#include <gsl/gsl_integration.h>

// The nine NLO impact-factor coefficient functions of arXiv:2310.06640:
// I2, J1, H2, K3 (real corrections, various qq/qg/gq/gg channels), J, Jv,
// Jv2 (Nc-channel real/virtual pair), H3, H5 (gg-channel real corrections).
// Each numerically integrates an underlying x,phi double integral with a
// channel-specific integrand selected by an internal flag (see
// nlo_coefficients.cpp). JJv_xi1 is the xi->1 subtraction term used only
// when params::with_xi1 is enabled. All depend on rp.alpha_s_running (the
// running-coupling prescription) and are evaluated against a fixed-Y
// dipole-amplitude slice (DipoleAmplitude1DSlice).
//
// Owns the GSL integration workspaces the double integral needs -- built
// once per instance and reused across every I2/J/.../JJv_xi1 call on it,
// rather than the file-scope-global/threadprivate pattern this replaced.
// Construct one per rapidity slice (matching how often the original code
// allocated its workspaces).
class NLOCoefficients{
public:
  explicit NLOCoefficients(const DipoleAmplitude& dipole);
  ~NLOCoefficients();
  NLOCoefficients(const NLOCoefficients&) = delete;
  NLOCoefficients& operator=(const NLOCoefficients&) = delete;

  double I2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double J(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double J1(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double H2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double H3(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double H5(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double K3(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double Jv(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double Jv2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  double JJv_xi1(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r);

private:
  double minlnr_, maxlnr_;
  gsl_integration_workspace *w_x_, *w_phi_;

  double func(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d,
              double r, double xi, double flag);

  static double integrand_x(double lnx, void *userdata);
  static double integrand_phi(double phi, void *userdata);
};
