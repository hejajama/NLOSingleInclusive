#pragma once

#include "Sr_interp_1D.hpp"
#include "dipole_amplitude.hpp"
#include "params.hpp"
#include "spline_wrappers.hpp"

#include <array>

class PdfSet;   // pdf_set.hpp -- only used by reference here

// Everything that depends on a single kinematic point (xp, xg, k): bicubic
// (r,Y) interpolations of the ten NLO coefficient functions
// (nlo_coefficients.hpp -- see that header for the paper-symbol mapping,
// which is not 1:1), their xi-convolution with the PDF at this point, and
// the fixed-Y dipole-amplitude slice used by the LO cross section
// (sigma_LO.cpp).
//
// xi_convolution(r) is the parton-level NLO correction to the cross
// section at fixed r: the dxi integral in arXiv:2310.06640 Eqs. (9)
// (qq channel) and (11a-11c) (gg/qg/gq channels), sec. II -- see
// integrand_xi() in point_tables.cpp for the channel-by-channel
// breakdown (params::RunParameters::channel selects which of Eq. (9)
// [Channel::QQ] / Eq. (11a) [Channel::GG] / Eq. (11b) [Channel::QG] /
// Eq. (11c) [Channel::GQ] applies).
//
// Built once per (z,k) point in main()'s loop; replaces the old
// init_interp()/clear_interp() free-function pair with a
// constructor/destructor.
class PointTables{
public:
  PointTables(const params::RunParameters& rp, const PdfSet& pdf,
              const DipoleAmplitude& dipole, double xp, double xg, double k);
  PointTables(const PointTables&) = delete;
  PointTables& operator=(const PointTables&) = delete;

  double I2(double r, double y) const;
  double J(double r, double y) const;
  double K1(double r, double y) const;
  double H2(double r, double y) const;
  double H3(double r, double y) const;
  double H4(double r, double y) const;
  double K2(double r, double y) const;
  double Jv(double r, double y) const;
  double Jv2(double r, double y) const;
  double JJv_xi1(double r, double y) const;

  // The xi-convolution of the coefficient functions with the PDF at this
  // point's (xp, xg, k).
  double xi_convolution(double r) const;

  // S(r, Y) at this point's fixed rapidity Y = log(1/xg).
  double dipole_slice(double r) const;

private:
  enum Coeff{ kI2, kJ, kK1, kH2, kH3, kH4, kK2, kJv, kJv2, kJJv_xi1, kNumCoeffs };

  double minr_, maxr_;   // for xi_convolution()'s bounds clamp
  double y_floor_;       // log(1/dipole.x0()): floor applied to y before
                         // evaluating any coefficient-function spline
  std::array<Spline2D, kNumCoeffs> coeff_splines_;
  Spline1D xi_conv_spline_;
  DipoleAmplitude1DSlice outer_slice_;

  void build_coefficient_tables(const params::RunParameters& rp,
                                 const DipoleAmplitude& dipole, double xg);
  void build_xi_convolution(const params::RunParameters& rp, const PdfSet& pdf,
                             const DipoleAmplitude& dipole, double xp, double xg, double k);

  static double integrand_xi(double xi, void *userdata);
};
