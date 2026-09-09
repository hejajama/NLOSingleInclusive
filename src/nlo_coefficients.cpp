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
  // (integrand_phi sets phi; integrand_x reads everything).
  struct IntegrandContext{
    NLOCoefficients *self;
    const RunParameters *rp;
    const DipoleAmplitude1DSlice *sr1d;
    double r;
    double xi;
    double phi;
    double flag;
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
  double flag=ctx.flag;
  // flag=1: J, flag=-1: Jv, flag=-2: Jv2, flag=0: J-Jv(xi=1), flag=2: I2,
  // flag=3: K1/4 (the K1 wrapper multiplies by 4, see below), flag=4: H2,
  // flag=5: H3, flag=6: H4, flag=7: K2/4 (K2 wrapper multiplies by 4)

  double x=exp(lnx), x2=Sq(x), r2=Sq(r), sprx=r*x*cos(phi), rpx2=r2+x2+2*sprx, rmx2=r2+x2-2*sprx;

  double res;
  if(flag>6.5){
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
  }
  else if(flag>5.5){
    // H4, Eq. (12d)
    double dip1 = sqrt(r2 + Sq(xi)*x2 - 2*xi*sprx);
    double dip2 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    res = 2 * (sr1d(dip2) * sr1d(dip1) - Sq(sr1d(dip2)));
    if(rp.alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
  }
  else if(flag>4.5){
    // H3
    double dip1 = sqrt(r2 + Sq(xi)*x2 + 2*xi*sprx);
    double dip2 = sqrt(r2 + Sq(1-xi)*x2 - 2*(1-xi)*sprx);
    double dip3 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    res = 4 * (sr1d(x) * sr1d(dip1) * sr1d(dip2) - Sq(sr1d(dip3)));
    if(rp.alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
  }
  else if(flag>3.5){
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
  }
  else if(flag>2.5){
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
  }
  else if(flag>1.5){
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

  }else if(flag>0.5){
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

  }else if(flag>-0.5){
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
  }else if(flag>-1.5){
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

  }else{
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

    }
  //if(alpha_s_running==DAUGHTER){
  //  res*=alpha_s_pos(x);
  //}
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
                              double r, double xi, double flag){
  IntegrandContext ctx{this, &rp, &sr1d, r, xi, 0.0, flag};
  double result, error;
  gsl_function F;
  F.function=&NLOCoefficients::integrand_phi;
  F.params=&ctx;
  gsl_integration_qag(&F,0,M_PI,epsabs_gsl,epsrel_gsl,gsl_maxpoints,
                      GSL_INTEG_GAUSS15,w_phi_,&result,&error);
  return (2*result)/Sq(2*M_PI); // 2: int over pi instead of 2*pi
}


double NLOCoefficients::I2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,2);
}


double NLOCoefficients::J(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 2*func(rp,sr1d,r,xi,1);
}


double NLOCoefficients::K1(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 4*func(rp,sr1d,r,xi,3);
}


double NLOCoefficients::H2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,4);
}


double NLOCoefficients::H3(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,5);
}


double NLOCoefficients::H4(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return func(rp,sr1d,r,xi,6);
}


double NLOCoefficients::K2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 4*func(rp,sr1d,r,xi,7);
}


double NLOCoefficients::Jv(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 2*func(rp,sr1d,r,xi,-1);
}


double NLOCoefficients::Jv2(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r, double xi){
  return 2*func(rp,sr1d,r,xi,-2);
}


double NLOCoefficients::JJv_xi1(const RunParameters& rp, const DipoleAmplitude1DSlice& sr1d, double r){
  return func(rp,sr1d,r,1,0);
}
