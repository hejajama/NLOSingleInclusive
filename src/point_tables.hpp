#pragma once

#include "Sr_interp_1D.hpp"
#include "dipole_amplitude.hpp"
#include "params.hpp"
#include "spline_wrappers.hpp"

#include <array>

class PdfSet;   // pdf_set.hpp -- only used by reference here

// Everything that depends on a single kinematic point (xp, xg, k): bicubic
// (r,Y) interpolations of the ten NLO coefficient functions
// (nlo_coefficients.hpp), their xi-convolution with the PDF at this point
// (arXiv:2310.06640 sec. 3), and the fixed-Y dipole-amplitude slice used by
// the LO cross section (sigma_LO.cpp).
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
  double J1(double r, double y) const;
  double H2(double r, double y) const;
  double H3(double r, double y) const;
  double H5(double r, double y) const;
  double K3(double r, double y) const;
  double Jv(double r, double y) const;
  double Jv2(double r, double y) const;
  double JJv_xi1(double r, double y) const;

  // The xi-convolution of the coefficient functions with the PDF at this
  // point's (xp, xg, k).
  double xi_convolution(double r) const;

  // S(r, Y) at this point's fixed rapidity Y = log(1/xg).
  double dipole_slice(double r) const;

private:
  enum Coeff{ kI2, kJ, kJ1, kH2, kH3, kH5, kK3, kJv, kJv2, kJJv_xi1, kNumCoeffs };

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
