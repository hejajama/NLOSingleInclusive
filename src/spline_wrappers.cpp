#include "spline_wrappers.hpp"

Spline1D::~Spline1D(){
  if(spline_) gsl_spline_free(spline_);
  if(acc_) gsl_interp_accel_free(acc_);
}

void Spline1D::build(const double *xvals, const double *yvals, int n){
  acc_ = gsl_interp_accel_alloc();
  spline_ = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline_, xvals, yvals, n);
}

double Spline1D::eval(double x) const{
  return gsl_spline_eval(spline_, x, acc_);
}


Spline2D::~Spline2D(){
  if(spline_) gsl_spline2d_free(spline_);
  if(xacc_) gsl_interp_accel_free(xacc_);
  if(yacc_) gsl_interp_accel_free(yacc_);
}

void Spline2D::init(int rpoints, int ypoints){
  rpoints_ = rpoints;
  ypoints_ = ypoints;
  backing_.resize(static_cast<std::size_t>(rpoints) * ypoints);
  spline_ = gsl_spline2d_alloc(gsl_interp2d_bicubic, rpoints, ypoints);
  xacc_ = gsl_interp_accel_alloc();
  yacc_ = gsl_interp_accel_alloc();
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
  return gsl_spline2d_eval(spline_, r, y, xacc_, yacc_);
}
