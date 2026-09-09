#include "nlo_coefficients.hpp"

#include "common.hpp"
#include "running_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <gsl/gsl_math.h>

using namespace std;
using namespace params;

namespace {
  // Context threaded through the phi/x double integral's GSL callbacks
  // (integrand_phi sets phi; integrand_x reads everything). `term`
  // selects the branch integrand_x() takes -- see NLOCoefficients::Term
  // in nlo_coefficients.hpp for what each one computes.
  struct IntegrandContext{
    NLOCoefficients *self;
    const RunParameters *rp;
    const DipoleAmplitude1DSlice *sr1d;
    double r;
    double xi;
    double phi;
    NLOCoefficients::Term term;
  };
}


NLOCoefficients::NLOCoefficients(const DipoleAmplitude& dipole)
  : minlnr_(dipole.min_ln_r()), maxlnr_(dipole.max_ln_r())
{
  w_x_ = gsl_integration_workspace_alloc(gsl_maxpoints);
  w_phi_ = gsl_integration_workspace_alloc(gsl_maxpoints);
}


NLOCoefficients::~NLOCoefficients(){
  gsl_integration_workspace_free(w_x_);
  gsl_integration_workspace_free(w_phi_);
}


double NLOCoefficients::integrand_x(double lnx, void *userdata){
  const IntegrandContext &ctx = *static_cast<IntegrandContext*>(userdata);
  const RunParameters &rp = *ctx.rp;
  const DipoleAmplitude1DSlice &sr1d = *ctx.sr1d;
  double r=ctx.r;
  double xi=ctx.xi;
  double phi=ctx.phi;
  Term term=ctx.term;

  double x=exp(lnx), x2=Sq(x), r2=Sq(r), sprx=r*x*cos(phi), rpx2=r2+x2+2*sprx, rmx2=r2+x2-2*sprx;

  double res;
  switch(term){
  case Term::K2: {
    // K2/4 -- the K2() wrapper multiplies this by 4 to match Eq. (12f) exactly
    double dip1 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    double K = (x2+sprx)/(x2*rpx2);
    res = 2 * exp(2*lnx) * K * sr1d(xi*x) * sr1d(dip1);
    if(rp.alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
    else if(rp.alpha_s_running==MIXED){
      res *= alpha_s_pos(min(xi*x, xi*sqrt(rpx2)));
    }
    else if(rp.alpha_s_running==PARENT){
      res *= alpha_s_pos(r);
    }
    break;
  }
  case Term::H4: {
    // H4, Eq. (12d)
    double dip1 = sqrt(r2 + Sq(xi)*x2 - 2*xi*sprx);
    double dip2 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    res = 2 * (sr1d(dip2) * sr1d(dip1) - Sq(sr1d(dip2)));
    if(rp.alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
    break;
  }
  case Term::H3: {
    // H3
    double dip1 = sqrt(r2 + Sq(xi)*x2 + 2*xi*sprx);
    double dip2 = sqrt(r2 + Sq(1-xi)*x2 - 2*(1-xi)*sprx);
    double dip3 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    res = 4 * (sr1d(x) * sr1d(dip1) * sr1d(dip2) - Sq(sr1d(dip3)));
    if(rp.alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
    break;
  }
  case Term::H2: {
    // H2
    double dip1 = sqrt(Sq(xi)*rpx2);
    double dip2 = sqrt(Sq(xi)*r2 + Sq(1-xi)*x2 - 2*xi*(1-xi)*sprx);
    double K = (x2+sprx)/(x2*rpx2);
    res = 8 * exp(2*lnx) * K * sr1d(x) * sr1d(dip1) * sr1d(dip2);
    if(rp.alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
    else if(rp.alpha_s_running==PARENT){
        res *= alpha_s_pos(r);
    }
    else if(rp.alpha_s_running==MIXED){
        res *= alpha_s_pos(min(xi*x, xi*sqrt(rpx2)));
    }
    break;
  }
  case Term::K1: {
    // K1/4 -- the K1() wrapper multiplies this by 4 to match Eq. (12e) exactly
    double dip1 = sqrt(Sq(xi)*r2+Sq(1-xi)*x2-2*xi*(1-xi)*sprx);
    double K = (x2+sprx)/(x2*rpx2);
    res = 2 * exp(2*lnx) * K * sr1d(x) * sr1d(dip1);
      /*
    if(alpha_s_running==MIXED){
      res*=alpha_s_pos(min(xi*x,xi*sqrt(rpx2)));
    }
    else if(alpha_s_running==MIXEDBD){
      res*=alpha_s_pos(x);
    }
      */
    break;
  }
  case Term::I2: {
    // I2
    double dip1=sqrt(Sq(xi)*r2+Sq(1-xi)*x2-2*xi*(1-xi)*sprx);
    double dip2=sqrt(rpx2);
    double K=(x2+sprx)/(x2*rpx2);
    res=exp(2*lnx)*K*sr1d(dip1);
    if(rp.alpha_s_running==PARENT){
      res*=alpha_s_pos(r);
    }
    else if(rp.alpha_s_running==DAUGHTER || rp.alpha_s_running==SMALLEST){
      res*=alpha_s_pos(dip2);
    }
    break;
  }
  case Term::J: {
    // J
    double dip1=sqrt(r2+Sq(1-xi)*x2+2*(1-xi)*sprx);
    double dip2=xi*x;
    double dip3=sqrt(rpx2);
    double K=(x2+sprx)/(x2*rpx2);
    if(rp.alpha_s_running==PARENT){
      res=alpha_s_pos(r)*exp(2*lnx)*K*(sr1d(dip1)-sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==DAUGHTER){
      res=alpha_s_pos(x)*exp(2*lnx)*K*(sr1d(dip1)-sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==SMALLEST){
      res=exp(2*lnx)*K*(alpha_s_pos(x)*sr1d(dip1)-alpha_s_pos(min(r, min(xi*x, xi*dip3)))*sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==MIXED){
      res=exp(2*lnx)*K*alpha_s_pos(min(xi*x,xi*dip3))*(sr1d(dip1)-sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==MIXEDBD){
      res=exp(2*lnx)*K*alpha_s_pos(x)*(sr1d(dip1)-sr1d(dip2)*sr1d(dip3));
    }
    else{
      res=exp(2*lnx)*K*(sr1d(dip1)-sr1d(dip2)*sr1d(dip3));
    }
    break;
  }
  case Term::JJv_xi1: {
    // J-Jv(xi=1)
    double K;
    if(x<sr1d.min_r()/10 || sqrt(rpx2)<sr1d.min_r()/10){
      K=0;
    }else{
      K=r2/(x2*rpx2+1e-20);
    }
    res=exp(2*lnx)*K*(sr1d(x)*sr1d(sqrt(rpx2))-sr1d(r));
    if(rp.alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
    break;
  }
  case Term::Jv: {
    // Jv, now include Jv2 stuffs already
    double dip1=sqrt(r2+Sq(1-xi)*x2-2*(1-xi)*sprx);
    double dip2=x;
    double dip3=sqrt(Sq(xi)*x2+r2+2*xi*sprx);
    double K=1/x2;
    if(rp.alpha_s_running==PARENT){
      // take 1
      res=exp(2*lnx)*K*(alpha_s_pos(dip1)*sr1d(dip1) - alpha_s_pos(r)*sr1d(dip2)*sr1d(dip3));
      // take 2
      //res=exp(2*lnx)*K*alpha_s_pos(r)*(sr1d(dip1) - sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==DAUGHTER){
      res=exp(2*lnx)*K*(alpha_s_pos(dip1)*sr1d(dip1) - alpha_s_pos(x)*sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==SMALLEST){
      // take 1
      res=exp(2*lnx)*K*(alpha_s_pos(dip1)*sr1d(dip1) - alpha_s_pos(min(r, min(xi*x, dip3)))*sr1d(dip2)*sr1d(dip3));
      // take 2
      //res=exp(2*lnx)*K*alpha_s_pos(min(r, min(xi*x, dip3)))*(sr1d(dip1) - sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==MIXED){
      res=exp(2*lnx)*K*alpha_s_pos(min(xi*x,xi*x+r))*(sr1d(dip1) - sr1d(dip2)*sr1d(dip3));
    }
    else if(rp.alpha_s_running==MIXEDBD){
      res=exp(2*lnx)*K*alpha_s_pos(x)*(sr1d(dip1) - sr1d(dip2)*sr1d(dip3));
    }
    else{
      res = exp(2*lnx)*K*(sr1d(dip1) - sr1d(dip2)*sr1d(dip3));
    }
    break;
  }
  case Term::Jv2: {
    // Jv2
    double dip1=sqrt(r2+Sq(1-xi)*x2-2*(1-xi)*sprx);
    double dip2=x;
    double dip3=sqrt(Sq(xi)*x2+r2+2*xi*sprx);
    double K=1/x2;
    if(rp.alpha_s_running==PARENT){
      // take 1
      res=-exp(2*lnx)*K*alpha_s_pos(r)*sr1d(dip2)*sr1d(dip3);
      // take 2
      //res=-exp(2*lnx)*K*alpha_s_pos(r)*sr1d(dip2)*sr1d(dip3);
    }
    else if(rp.alpha_s_running==DAUGHTER){
      res=-exp(2*lnx)*K*alpha_s_pos(x)*sr1d(dip2)*sr1d(dip3);
    }
    else if(rp.alpha_s_running==SMALLEST){
      // take 1
      res=-exp(2*lnx)*K*alpha_s_pos(min(r, min(xi*x, dip3)))*sr1d(dip2)*sr1d(dip3);
      // take 2
      //res=-exp(2*lnx)*K*alpha_s_pos(min(r, min(xi*x, dip3)))*sr1d(dip2)*sr1d(dip3);
    }
    else if(rp.alpha_s_running==MIXED){
      res=-exp(2*lnx)*K*alpha_s_pos(min(xi*x, dip3))*sr1d(dip2)*sr1d(dip3);
    }
    else{
        res=-exp(2*lnx)*K*sr1d(dip2)*sr1d(dip3);
    }
    break;
  }
  }
  if(gsl_finite(res)==1){
    return res;
  }else{
    return 0;
  }
}


double NLOCoefficients::integrand_phi(double phi, void *userdata){
  IntegrandContext &ctx = *static_cast<IntegrandContext*>(userdata);
  ctx.phi=phi;
  NLOCoefficients &self = *ctx.self;
  double result, error;
  gsl_function F;
  F.function=&NLOCoefficients::integrand_x;
  F.params=userdata;
  gsl_integration_qag(&F,self.minlnr_,self.maxlnr_,epsabs_gsl,epsrel_gsl,gsl_maxpoints,
                      GSL_INTEG_GAUSS15,self.w_x_,&result,&error);
  return result;
}


double NLOCoefficients::func(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d,
                              double r, double xi, Term term){
  IntegrandContext ctx{this, &rp, &sr1d, r, xi, 0.0, term};
  double result, error;
  gsl_function F;
  F.function=&NLOCoefficients::integrand_phi;
  F.params=&ctx;
  gsl_integration_qag(&F,0,M_PI,epsabs_gsl,epsrel_gsl,gsl_maxpoints,
                      GSL_INTEG_GAUSS15,w_phi_,&result,&error);
  return (2*result)/Sq(2*M_PI); // 2: int over pi instead of 2*pi
}


double NLOCoefficients::I2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,Term::I2);
}


double NLOCoefficients::J(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 2*func(rp,sr1d,r,xi,Term::J);
}


double NLOCoefficients::K1(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 4*func(rp,sr1d,r,xi,Term::K1);
}


double NLOCoefficients::H2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,Term::H2);
}


double NLOCoefficients::H3(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,Term::H3);
}


double NLOCoefficients::H4(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,Term::H4);
}


double NLOCoefficients::K2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 4*func(rp,sr1d,r,xi,Term::K2);
}


double NLOCoefficients::Jv(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 2*func(rp,sr1d,r,xi,Term::Jv);
}


double NLOCoefficients::Jv2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 2*func(rp,sr1d,r,xi,Term::Jv2);
}


double NLOCoefficients::JJv_xi1(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r){
  return func(rp,sr1d,r,1,Term::JJv_xi1);
}
