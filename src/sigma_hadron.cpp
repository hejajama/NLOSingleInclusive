// Choice of n_zpoints's default (16, cli.hpp's --z-points): checked
// against pp, incoming=q, outgoing=q, rc=mom, p_h=2, zmin=0.4 (--col pp
// --b 0 --incoming q --outgoing q --rc mom --p 2 --muratio 1 --level
// hadron --zmin 0.4). LO (cheap, no PointTables) converges essentially
// immediately: n=2 is within 0.05% of an n=128 reference. NLO is pricier
// to check (each point costs a full PointTables build) but doubling
// n=8->16 moved it by 0.76% (0.000177201 -> 0.000175871) -- right at the
// edge of a 1% tolerance -- so 16 (rather than 8) is used as the default,
// for margin against whatever residual error remains at that point (the
// n=8->16 step already suggests it's well under 1%). Re-check this if the
// FF set, channel, or kinematic range changes enough that the z-integrand
// could develop sharper features.
#include "sigma_hadron.hpp"

#include "common.hpp"
#include "point_tables.hpp"
#include "sigma_LO.hpp"
#include "sigma_NLO.hpp"

#include <cmath>
#include <gsl/gsl_integration.h>

using namespace params;

HadronSigma sigma_hadron_ph(const RunParameters& rp, const PdfSet& pdf, const FfSet& ff,
                             const DipoleAmplitude& dipole, double p_h, double zmin,
                             double sqrts, double y, int n_zpoints){
  gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(n_zpoints);

  HadronSigma result{0, 0};
  for(int i=0; i<n_zpoints; i++){
    double z, w;
    gsl_integration_glfixed_point(zmin, 1, i, &z, &w, table);

    double k = p_h/z;
    double xp = (k/sqrts)*exp(y);    // Eq. 6, at this z's k
    double xg = (k/sqrts)*exp(-y);   // Eq. 5, at this z's k
    double D = ff.zD(rp, z, rp.mu2_ff);
    double jacobian = w/Sq(z);

    // LO only ever needs S(r,Y) (sigma_LO.cpp), never PointTables' NLO
    // coefficient/xi-convolution tables, so it skips PointTables entirely.
    result.LO += jacobian*D*sigma_LO_k(rp, pdf, dipole, xg, k, xp);

    // NLO does need the full PointTables, built once here and reused for
    // nothing else -- each Gauss-Legendre node is independent by
    // construction, so there's nothing to cache across them.
    PointTables tables(rp, pdf, dipole, xp, xg, k);
    result.NLO += jacobian*D*sigma_NLO_k(rp, pdf, tables, k, xp);
  }

  gsl_integration_glfixed_table_free(table);
  return result;
}
