#include "sigma_LO.hpp"

#include "common.hpp"
#include "intde1.hpp"

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

using namespace params;

namespace {
  // Eq. 7a (quark channel)/7b (gluon channel), given S(r,Y) (or S(r,Y)^2)
  // at this point.
  double sigma_LO_r_from_S(const RunParameters& rp, const PdfSet& pdf, double xp, double Sr){
    if(rp.channel==Channel::QQ){
      return pdf.xf(rp,xp,rp.mu2)*Sr;        // Eq. 7a, quark channel
    }
    else if(rp.channel==Channel::GG){
      return pdf.xf(rp,xp,rp.mu2)*Sq(Sr);    // Eq. 7b, gluon channel
    }
    return 0;
  }

  // Context threaded through integrand_sigma_LO_k via intdeo's userdata
  // pointer (intdeo's f has no userdata slot of its own to smuggle k/xp
  // through otherwise, unlike GSL's own gsl_function convention).
  struct SigmaLOContext{
    const RunParameters *rp;
    const PdfSet *pdf;
    const DipoleAmplitude *dipole;
    double Y;
    double k;
    double xp;
  };

  double integrand_sigma_LO_k(double r, void *userdata){
    const SigmaLOContext &ctx = *static_cast<SigmaLOContext*>(userdata);
    double Sr = ctx.dipole->S(r, ctx.Y);
    return r*gsl_sf_bessel_J0(ctx.k*r)*sigma_LO_r_from_S(*ctx.rp,*ctx.pdf,ctx.xp,Sr);
  }
}


double sigma_LO_k(const RunParameters& rp, const PdfSet& pdf, const DipoleAmplitude& dipole,
                   double xg, double k, double xp){
  SigmaLOContext ctx{&rp, &pdf, &dipole, std::log(1/xg), k, xp};
  double integral, error;
  intdeo(integrand_sigma_LO_k,0,k,epsrel_intde,&integral,&error,&ctx);
  return integral/(2*M_PI);
}
