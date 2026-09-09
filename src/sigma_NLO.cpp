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
  if(rp.channel==Channel::QQ){
      res += pdf.xf(rp,xp,rp.mu2)*Sr_0(rp,r);
  }
  if(rp.channel==Channel::GG){
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
  // rp.channel is always exactly one of QQ/QG/GQ/GG, so this condition is a
  // tautology (with_xi1's dead-code case never has to carry it alone) --
  // kept in this form, matching the original with_Nc/with_CF/with_gl/
  // with_gq/with_gg OR chain, rather than collapsed to `true`.
  if(rp.channel==Channel::QQ || with_xi1 || rp.channel==Channel::QG ||
     rp.channel==Channel::GQ || rp.channel==Channel::GG) res+=as*tables.xi_convolution(r);
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
