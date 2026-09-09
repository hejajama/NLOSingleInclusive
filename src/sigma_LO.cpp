#include "sigma_LO.hpp"

#include "Sr_interp_1D.hpp"
#include "common.hpp"
#include "intde1.hpp"
#include "params.hpp"
#include "utils.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

using namespace params;

namespace {
  // Context threaded through integrand_sigma_LO_k via intdeo's userdata
  // pointer (intdeo's f has no userdata slot of its own to smuggle k/xp
  // through otherwise, unlike GSL's own gsl_function convention).
  struct SigmaLOContext{
    const RunParameters *rp;
    double k;
    double xp;
  };
}


double sigma_LO_r(const RunParameters& rp, double r, double xp){
  if(rp.with_Nc||rp.with_CF){
    return xf(rp,xp,rp.mu2)*Sr_interp_1D(r);
  }
  else if(rp.with_gg){
    return xf(rp,xp,rp.mu2)*Sq(Sr_interp_1D(r));
  }
  return 0;
}


double integrand_sigma_LO_k(double r, void *userdata){
  const SigmaLOContext &ctx = *static_cast<SigmaLOContext*>(userdata);
  return r*gsl_sf_bessel_J0(ctx.k*r)*sigma_LO_r(*ctx.rp,r,ctx.xp);
}


double sigma_LO_k(const RunParameters& rp, double k, double xp){
  SigmaLOContext ctx{&rp, k, xp};
  double integral, error;
  intdeo(integrand_sigma_LO_k,0,k,epsrel_intde,&integral,&error,&ctx);
  return integral/(2*M_PI);
}

/*
double integrand_sigma_LO_p(double z){
    using namespace params;
    using namespace sigma_LO_p_tmp;

    // Compute FF(z)
    fragini_.fini=0;
    double FF;
    int ih_tmp = ih;
    int ic_tmp = ic;
    int io_tmp = io;
    int parton_tmp = parton;
    double scalesqr_tmp = scalesqr;
    fdss_(ih_tmp, ic_tmp, io_tmp, z, scalesqr_tmp, parton_tmp, FF);

    // Compute sigma_LO_k(p/z, tau/z)
    double k = p_tmp/z;
    mu2=Sq(k);
    double xp=(k/SQRTS)*exp(yh);
    double xg=(k/SQRTS)*exp(-yh);
    init_interp(xp,xg,k);
    double mult = sigma_LO_k(k,xp);
    clear_interp();

    return mult*FF/Sq(z);
}


double sigma_LO_p(double p){
    using namespace sigma_LO_p_tmp;
    p_tmp = p;
    double integral, error;
    double zmin = 0.05;
    intde(integrand_sigma_LO_p, zmin, 1, epsrel_intde, &integral, &error);
    return integral;
}
*/
