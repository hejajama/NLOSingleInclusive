#pragma once

#include "Sr_interp_1D.hpp"
#include "dipole_amplitude.hpp"
#include "params.hpp"

#include <gsl/gsl_integration.h>

// The NLO impact-factor coefficient functions of arXiv:2310.06640 sec. II
// that require a two-dimensional transverse integral (over the emitted
// parton's momentum x_perp and the angle phi between x_perp and r_perp) --
// see nlo_coefficients.cpp for the shared x,phi double integral and the
// channel-specific integrand each one selects via the Term enum below.
// (I1 and H1, Eqs. (10c)/(12a), need no such integral -- they're plain
// algebraic functions of S(r,Y), and live as free functions in
// point_tables.cpp instead.)
//
// Named to match the paper's own symbols (K1/K2/H4) -- these used to be
// called J1/K3/H5, an unrelated earlier naming with no connection to the
// paper's own K1/K2/H4 numbering (K1/K2 also needed rescaling by 4 to
// match the paper's literal definitions, not just renaming; see
// docs/PAPER_MAPPING.md, "Former J1/K3/H5 naming", for the history and
// derivation). See docs/PAPER_MAPPING.md for how each one plugs into the
// qq/qg/gq/gg channel cross sections (Eqs. 9, 11).
//
// Jv2 and JJv_xi1 are not part of the "unsubtracted scheme" this paper
// uses (Eqs. 9-12); JJv_xi1 is gated behind params::with_xi1, which is
// permanently false in this codebase, so it is never actually evaluated.
// Both appear to be leftover machinery for the alternative "subtracted
// scheme" of Ducloué et al., Phys. Rev. D 97 (2018) 054020
// [arXiv:1712.07480] (ref. [26] in arXiv:2310.06640) rather than anything
// defined in this paper -- unverified, since that paper's equations
// weren't checked for this pass.
//
// Owns the GSL integration workspaces the double integral needs -- built
// once per instance and reused across every I2/J/.../JJv_xi1 call on it,
// rather than the file-scope-global/threadprivate pattern this replaced.
// Construct one per rapidity slice (matching how often the original code
// allocated its workspaces).
class NLOCoefficients{
public:
  // Which term func()'s shared phi/x double integral evaluates -- see
  // nlo_coefficients.cpp for what each one computes. Exposed here only
  // because func()'s declaration below needs the type; nothing outside
  // this class and nlo_coefficients.cpp should need to name it.
  enum class Term{ Jv2 = -2, Jv = -1, JJv_xi1 = 0, J = 1, I2 = 2, K1 = 3, H2 = 4, H3 = 5, H4 = 6, K2 = 7 };

  explicit NLOCoefficients(const DipoleAmplitude& dipole);
  ~NLOCoefficients();
  NLOCoefficients(const NLOCoefficients&) = delete;
  NLOCoefficients& operator=(const NLOCoefficients&) = delete;

  // qq-channel real correction, Eq. (10d).
  double I2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // qq-channel real correction, Eq. (10a).
  double J(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // qg-channel term, Eq. (12e).
  double K1(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gg-channel term, Eq. (12b) -- same symbol as the paper.
  double H2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gg-channel term, Eq. (12c) -- same symbol as the paper.
  double H3(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gg-channel term, Eq. (12d).
  double H4(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gq-channel term, Eq. (12f).
  double K2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // qq-channel virtual correction, Eq. (10b).
  double Jv(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // Not part of Eqs. (9)-(12); see the class-level comment above.
  double Jv2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // Gated behind params::with_xi1 (always false); see the class-level
  // comment above.
  double JJv_xi1(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r);

private:
  double minlnr_, maxlnr_;
  gsl_integration_workspace *w_x_, *w_phi_;

  double func(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d,
              double r, double xi, Term term);

  static double integrand_x(double lnx, void *userdata);
  static double integrand_phi(double phi, void *userdata);
};
