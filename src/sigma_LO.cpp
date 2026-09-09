#include "sigma_LO.hpp"

#include "common.hpp"
#include "intde1.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

using namespace params;

namespace {
  // Context threaded through integrand_sigma_LO_k via intdeo's userdata
  // pointer (intdeo's f has no userdata slot of its own to smuggle k/xp
  // through otherwise, unlike GSL's own gsl_function convention).
  struct SigmaLOContext{
    const RunParameters *rp;
    const PdfSet *pdf;
    const PointTables *tables;
    double k;
    double xp;
  };
}


double sigma_LO_r(const RunParameters& rp, const PdfSet& pdf, const PointTables& tables,
                   double r, double xp){
  if(rp.channel==Channel::QQ){
    return pdf.xf(rp,xp,rp.mu2)*tables.dipole_slice(r);   // Eq. 7a, quark channel
  }
  else if(rp.channel==Channel::GG){
    return pdf.xf(rp,xp,rp.mu2)*Sq(tables.dipole_slice(r));   // Eq. 7b, gluon channel
  }
  return 0;
}


double integrand_sigma_LO_k(double r, void *userdata){
  const SigmaLOContext &ctx = *static_cast<SigmaLOContext*>(userdata);
  return r*gsl_sf_bessel_J0(ctx.k*r)*sigma_LO_r(*ctx.rp,*ctx.pdf,*ctx.tables,r,ctx.xp);
}


double sigma_LO_k(const RunParameters& rp, const PdfSet& pdf, const PointTables& tables,
                   double k, double xp){
  SigmaLOContext ctx{&rp, &pdf, &tables, k, xp};
  double integral, error;
  intdeo(integrand_sigma_LO_k,0,k,epsrel_intde,&integral,&error,&ctx);
  return integral/(2*M_PI);
}
