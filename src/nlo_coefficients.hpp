#pragma once

#include "Sr_interp_1D.hpp"
#include "dipole_amplitude.hpp"
#include "params.hpp"

#include <gsl/gsl_integration.h>

// The NLO impact-factor coefficient functions of arXiv:2310.06640 sec. II
// that require a two-dimensional transverse integral (over the emitted
// parton's momentum x_perp and the angle phi between x_perp and r_perp) --
// see nlo_coefficients.cpp for the shared x,phi double integral and the
// channel-specific integrand each one selects via an internal flag.
// (I1 and H1, Eqs. (10c)/(12a), need no such integral -- they're plain
// algebraic functions of S(r,Y), and live as free functions in
// point_tables.cpp instead.)
//
// IMPORTANT: several of these functions are named differently here than
// in the paper -- the mapping is (paper symbol -> code symbol):
//   J   (Eq. 10a) -> J          Jv  (Eq. 10b) -> Jv         I2 (Eq. 10d) -> I2
//   H1  (Eq. 12a) -> H1 (point_tables.cpp)     H2 (Eq. 12b) -> H2
//   H3  (Eq. 12c) -> H3          H4  (Eq. 12d) -> H5   (NOT H4!)
//   K1  (Eq. 12e) -> 4*J1        K2  (Eq. 12f) -> 4*K3 (NOT K1/K2!)
// i.e. J1(r,xi) == (1/4) K1_paper(r,xi) and K3(r,xi) == (1/4) K2_paper(r,xi)
// -- the code's own "J1"/"K3" bear no relation to the paper's K1/K2
// numbering. See docs/PAPER_MAPPING.md for the full derivation and how
// each one plugs into the qq/qg/gq/gg channel cross sections (Eqs. 9, 11).
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
  explicit NLOCoefficients(const DipoleAmplitude& dipole);
  ~NLOCoefficients();
  NLOCoefficients(const NLOCoefficients&) = delete;
  NLOCoefficients& operator=(const NLOCoefficients&) = delete;

  // qq-channel real correction, Eq. (10d).
  double I2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // qq-channel real correction, Eq. (10a).
  double J(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // qg-channel term: J1 == (1/4) x paper's K1, Eq. (12e).
  double J1(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gg-channel term, Eq. (12b) -- same symbol as the paper.
  double H2(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gg-channel term, Eq. (12c) -- same symbol as the paper.
  double H3(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gg-channel term: H5 == paper's H4, Eq. (12d).
  double H5(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
  // gq-channel term: K3 == (1/4) x paper's K2, Eq. (12f).
  double K3(const params::RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi);
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
              double r, double xi, double flag);

  static double integrand_x(double lnx, void *userdata);
  static double integrand_phi(double phi, void *userdata);
};
