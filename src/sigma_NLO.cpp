#include "sigma_NLO.hpp"

#include "bksol_nlodisfit.hpp"
#include "common.hpp"
#include "intde1.hpp"
#include "params.hpp"
#include "utils.hpp"
#include "xi_int.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

using namespace params;

namespace {
  // Context threaded through integrand_sigma_NLO_k via intdeo's userdata
  // pointer (see sigma_LO.cpp for why this is needed at all).
  struct SigmaNLOContext{
    const RunParameters *rp;
    double k;
    double xp;
  };
}


double sigma_NLO_r(const RunParameters& rp, double r, double k, double xp){
  double res = 0;
  if(rp.with_Nc || rp.with_CF){
      res += xf(rp,xp,rp.mu2)*Sr_0(rp,r);
  }
  if(rp.with_gg){
      res += xf(rp,xp,rp.mu2)*Sq(Sr_0(rp,r));
  }
  double as;
  if(rp.alpha_s_running==FIXED){
    as=alpha_s_fixed;
  }else if(rp.alpha_s_running==MOM){
    as=alpha_s_mom(rp.mu2);
  }else{
    as=1;
  }
   if(rp.with_Nc || with_xi1 || rp.with_CF || rp.with_gl || rp.with_gq || rp.with_gg) res+=as*xi_int_interp(r);
  return res;
}


double integrand_sigma_NLO_k(double r, void *userdata){
  const SigmaNLOContext &ctx = *static_cast<SigmaNLOContext*>(userdata);
  return r*gsl_sf_bessel_J0(ctx.k*r)*sigma_NLO_r(*ctx.rp,r,ctx.k,ctx.xp);
}


double sigma_NLO_k(const RunParameters& rp, double k, double xp){
  SigmaNLOContext ctx{&rp, k, xp};
  double integral, error;
  intdeo(integrand_sigma_NLO_k,0,k,epsrel_intde,&integral,&error,&ctx);
  return integral/(2*M_PI);
}

/*
double integrand_sigma_NLO_p(double z){
    using namespace params;
    using namespace sigma_NLO_p_tmp;

    // Compute FF(z)
    fragini_.fini=0;
    double FF;
    int ih_tmp = ih;
    int ic_tmp = ic;
    int io_tmp = io;
    int parton_tmp = parton;
    double scalesqr_tmp = scalesqr;
    fdss_(ih_tmp, ic_tmp, io_tmp, z, scalesqr_tmp, parton_tmp, FF);

    // Compute sigma_NLO_k(p/z, tau/z)
    double k = p_tmp/z;
    mu2=Sq(k);
    double xp=(k/SQRTS)*exp(yh);
    double xg=(k/SQRTS)*exp(-yh);
    init_interp(xp,xg,k);
    double mult = sigma_NLO_k(k,xp);
    clear_interp();

    return mult*FF/Sq(z);
}


double sigma_NLO_p(double p){
    using namespace sigma_NLO_p_tmp;
    p_tmp = p;
    double integral, error;
    double zmin = 0.05;
    intde(integrand_sigma_NLO_p, zmin, 1, epsrel_intde, &integral, &error);
    return integral;
}
*/
