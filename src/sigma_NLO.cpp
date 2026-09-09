#include "sigma_NLO.hpp"

#include "common.hpp"
#include "dipole_amplitude.hpp"
#include "intde1.hpp"
#include "running_coupling.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

using namespace params;

namespace {
  // Context threaded through integrand_sigma_NLO_k via intdeo's userdata
  // pointer (see sigma_LO.cpp for why this is needed at all).
  struct SigmaNLOContext{
    const RunParameters *rp;
    const PdfSet *pdf;
    const PointTables *tables;
    double k;
    double xp;
  };
}


double sigma_NLO_r(const RunParameters& rp, const PdfSet& pdf, const PointTables& tables,
                    double r, double k, double xp){
  double res = 0;
  if(rp.with_Nc || rp.with_CF){
      res += pdf.xf(rp,xp,rp.mu2)*Sr_0(rp,r);
  }
  if(rp.with_gg){
      res += pdf.xf(rp,xp,rp.mu2)*Sq(Sr_0(rp,r));
  }
  double as;
  if(rp.alpha_s_running==FIXED){
    as=alpha_s_fixed;
  }else if(rp.alpha_s_running==MOM){
    as=alpha_s_mom(rp.mu2);
  }else{
    as=1;
  }
   if(rp.with_Nc || with_xi1 || rp.with_CF || rp.with_gl || rp.with_gq || rp.with_gg) res+=as*tables.xi_convolution(r);
  return res;
}


double integrand_sigma_NLO_k(double r, void *userdata){
  const SigmaNLOContext &ctx = *static_cast<SigmaNLOContext*>(userdata);
  return r*gsl_sf_bessel_J0(ctx.k*r)*sigma_NLO_r(*ctx.rp,*ctx.pdf,*ctx.tables,r,ctx.k,ctx.xp);
}


double sigma_NLO_k(const RunParameters& rp, const PdfSet& pdf, const PointTables& tables,
                    double k, double xp){
  SigmaNLOContext ctx{&rp, &pdf, &tables, k, xp};
  double integral, error;
  intdeo(integrand_sigma_NLO_k,0,k,epsrel_intde,&integral,&error,&ctx);
  return integral/(2*M_PI);
}
