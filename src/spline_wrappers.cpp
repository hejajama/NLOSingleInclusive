#include "spline_wrappers.hpp"

namespace{
// One accelerator per thread per axis, not per spline instance: GSL
// accelerators just cache a last-looked-up index and aren't tied to any
// particular spline's data, but mutating one from multiple threads at once
// is a race. thread_local gives each thread its own, freed at thread exit.
struct AccelGuard{
  gsl_interp_accel *p = gsl_interp_accel_alloc();
  ~AccelGuard(){ gsl_interp_accel_free(p); }
};
}

Spline1D::~Spline1D(){
  if(spline_) gsl_spline_free(spline_);
}

void Spline1D::build(const double *xvals, const double *yvals, int n){
  spline_ = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline_, xvals, yvals, n);
}

double Spline1D::eval(double x) const{
  thread_local AccelGuard acc;
  return gsl_spline_eval(spline_, x, acc.p);
}


Spline2D::~Spline2D(){
  if(spline_) gsl_spline2d_free(spline_);
}

void Spline2D::init(int rpoints, int ypoints){
  rpoints_ = rpoints;
  ypoints_ = ypoints;
  backing_.resize(static_cast<std::size_t>(rpoints) * ypoints);
  spline_ = gsl_spline2d_alloc(gsl_interp2d_bicubic, rpoints, ypoints);
}

void Spline2D::set(int ir, int iy, double value){
  gsl_spline2d_set(spline_, backing_.data(), ir, iy, value);
}

void Spline2D::build(const double *rvals, const double *yvals){
  gsl_spline2d_init(spline_, rvals, yvals, backing_.data(), rpoints_, ypoints_);
  backing_.clear();
  backing_.shrink_to_fit();
}

double Spline2D::eval(double r, double y) const{
  thread_local AccelGuard xacc, yacc;
  return gsl_spline2d_eval(spline_, r, y, xacc.p, yacc.p);
}
