#pragma once

#include <vector>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

// Small RAII wrappers around GSL's spline types, used throughout the
// dipole-amplitude and NLO-coefficient interpolation code instead of
// hand-paired alloc/free calls. Generic numerics helpers, not tied to any
// particular physics quantity.
//
// eval() is safe to call concurrently on the same instance from multiple
// threads: the underlying gsl_spline/gsl_spline2d data is immutable once
// build()/init()+build() have run, and eval() itself uses a thread_local
// gsl_interp_accel (spline_wrappers.cpp) rather than a per-instance one, so
// no accelerator is ever mutated by more than one thread.

// A 1D cubic spline built from a full (x,y) array in one shot.
class Spline1D{
public:
  Spline1D() = default;
  ~Spline1D();
  Spline1D(const Spline1D&) = delete;
  Spline1D& operator=(const Spline1D&) = delete;

  // xvals/yvals must have `n` entries each; xvals must be strictly
  // increasing (GSL's requirement for gsl_spline_init).
  void build(const double *xvals, const double *yvals, int n);
  double eval(double x) const;

private:
  gsl_spline *spline_ = nullptr;
};

// A bicubic 2D spline, filled in one value at a time (set()) and then
// finalized (build()) once the whole (rpoints x ypoints) grid is known --
// mirrors gsl_spline2d's own alloc/set/init/free lifecycle.
class Spline2D{
public:
  Spline2D() = default;
  ~Spline2D();
  Spline2D(const Spline2D&) = delete;
  Spline2D& operator=(const Spline2D&) = delete;

  void init(int rpoints, int ypoints);
  void set(int ir, int iy, double value);
  // rvals/yvals must have rpoints/ypoints entries respectively (as passed
  // to init()) and must be strictly increasing.
  void build(const double *rvals, const double *yvals);
  double eval(double r, double y) const;

private:
  gsl_spline2d *spline_ = nullptr;
  std::vector<double> backing_;
  int rpoints_ = 0, ypoints_ = 0;
};
