#pragma once

#include "dipole_amplitude.hpp"
#include "spline_wrappers.hpp"

// Cubic spline in r of S(r,Y) at a single, fixed rapidity Y -- a cheaper
// slice of DipoleAmplitude::S used inside r-integrals that run at fixed Y
// (the NLO coefficient functions in nlo_coefficients.cpp, and the LO cross
// section in sigma_LO.cpp).
class DipoleAmplitude1DSlice{
public:
  DipoleAmplitude1DSlice(const DipoleAmplitude& dipole, double Y);
  DipoleAmplitude1DSlice(const DipoleAmplitude1DSlice&) = delete;
  DipoleAmplitude1DSlice& operator=(const DipoleAmplitude1DSlice&) = delete;

  double operator()(double r) const;

  double min_r() const { return minr_; }

private:
  double minr_, maxr_;
  Spline1D spline_;
};
