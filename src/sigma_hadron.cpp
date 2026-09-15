#include "sigma_hadron.hpp"

#include "common.hpp"
#include "intde1.hpp"
#include "point_tables.hpp"
#include "sigma_LO.hpp"
#include "sigma_NLO.hpp"

#include <cmath>

using namespace params;

namespace {
  // Context threaded through the two integrand_*_ph functions via intde's
  // userdata pointer (see sigma_LO.cpp/sigma_NLO.cpp for why this pattern
  // is needed at all).
  struct SigmaHadronContext{
    const RunParameters *rp;
    const PdfSet *pdf;
    const FfSet *ff;
    const DipoleAmplitude *dipole;
    double p_h;
    double sqrts;
    double y;
  };

  // LO only ever needs S(r,Y) (sigma_LO.cpp's sigma_LO_r_from_S), never
  // PointTables' NLO coefficient/xi-convolution tables -- building a full
  // PointTables per z-node here would pay for machinery this integrand
  // never touches (point_tables.cpp's per-r xi-convolution and per-(r,y)
  // NLO-coefficient tables are what make that expensive), so it calls
  // dipole.S(r,Y) directly instead (sigma_LO.hpp).
  double integrand_sigma_LO_ph(double z, void *userdata){
    SigmaHadronContext &ctx = *static_cast<SigmaHadronContext*>(userdata);
    const RunParameters &rp = *ctx.rp;
    double k = ctx.p_h/z;
    double xp = (k/ctx.sqrts)*exp(ctx.y);    // Eq. 6, at this z's k
    double xg = (k/ctx.sqrts)*exp(-ctx.y);   // Eq. 5, at this z's k
    double D = ctx.ff->zD(rp, z, rp.mu2_ff);
    return D*sigma_LO_k(rp, *ctx.pdf, *ctx.dipole, xg, k, xp)/Sq(z);
  }

  // NLO does need the full PointTables (its coefficient/xi-convolution
  // tables), built fresh at every z-node: intde's own DE-quadrature
  // (intde1.cpp) samples an entirely new, disjoint set of transformed
  // node locations at each refinement level, so within this one z-integral
  // there's no repeated z to cache against.
  double integrand_sigma_NLO_ph(double z, void *userdata){
    SigmaHadronContext &ctx = *static_cast<SigmaHadronContext*>(userdata);
    const RunParameters &rp = *ctx.rp;
    double k = ctx.p_h/z;
    double xp = (k/ctx.sqrts)*exp(ctx.y);
    double xg = (k/ctx.sqrts)*exp(-ctx.y);
    PointTables tables(rp, *ctx.pdf, *ctx.dipole, xp, xg, k);
    double D = ctx.ff->zD(rp, z, rp.mu2_ff);
    return D*sigma_NLO_k(rp, *ctx.pdf, tables, k, xp)/Sq(z);
  }
}


HadronSigma sigma_hadron_ph(const RunParameters& rp, const PdfSet& pdf, const FfSet& ff,
                             const DipoleAmplitude& dipole, double p_h, double zmin,
                             double sqrts, double y){
  SigmaHadronContext ctx{&rp, &pdf, &ff, &dipole, p_h, sqrts, y};

  HadronSigma result;
  double error;
  intde(integrand_sigma_LO_ph, zmin, 1, epsrel_intde, &result.LO, &error, &ctx);
  intde(integrand_sigma_NLO_ph, zmin, 1, epsrel_intde, &result.NLO, &error, &ctx);
  return result;
}
