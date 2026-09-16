// Choice of n_zpoints's default (16, cli.hpp's --z-points): checked
// against pp, incoming=q, outgoing=q, rc=mom, p_h=2, zmin=0.4 (--col pp
// --b 0 --incoming q --outgoing q --rc mom --pt 2 --muratio 1 --level
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
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>
#include <gsl/gsl_integration.h>
#include <omp.h>

using namespace params;

namespace {
  // True for any rp.incoming that sums more than one flavor's PDF
  // together: the dedicated "q" (all Nf flavors) aggregate, or a '+'-joined
  // multi-flavor token like "u+ubar" (PdfSet::xf) -- as opposed to a single
  // explicit flavor ("u", "g", ...), which never has this problem.
  bool is_multiflavor_incoming(const std::string& incoming){
    return incoming=="q" || incoming.find('+')!=std::string::npos;
  }
}

HadronSigma sigma_hadron_ph(const RunParameters& rp, const PdfSet& pdf, const FfSet& ff,
                             const DipoleAmplitude& dipole, double p_h, double zmin,
                             double sqrts, double y, int n_zpoints){
  // rp.incoming=="q" (or a '+'-joined multi-flavor token, see
  // is_multiflavor_incoming above) makes PdfSet::xf sum the PDF over more
  // than one flavor; rp.outgoing=="q" independently makes FfSet::zD sum
  // the FF over all Nf flavors. The QQ channel's parton-level hard function
  // (sigma_LO_r_from_S/sigma_NLO_r, point_tables.cpp's QQ branches)
  // conserves flavor -- an incoming quark of flavor i stays flavor i
  // through to the fragmenting parton -- so the physically correct
  // hadron-level sum is sum_i f_i(x)*D_i(z), one term per matching flavor
  // pair. Multiplying a multi-flavor incoming sum by rp.outgoing=="q"'s
  // all-flavor FF sum instead computes (sum_i f_i)*(sum_j D_j): that
  // includes unphysical cross-flavor terms (e.g. an incoming u-quark
  // weighted by the d-quark FF) that don't correspond to any real process
  // -- UNLESS every flavor being summed on the incoming side happens to
  // share the exact same FF value, which rp.outgoing=="q" can't express
  // (it always means "all Nf flavors", not "whichever ones happen to
  // match"). So this combination is rejected outright rather than silently
  // returning a wrong number -- see README.md, "Hadron-level flavor sums",
  // for how to get the correct flavor-summed qq result (run once per
  // explicit matching flavor -- or flavor pair, e.g. "u+ubar" against a
  // single explicit outgoing "u" when the FF is known to match, see "pi0
  // fragmentation functions" -- and add the results).
  //
  // QG (incoming sums flavor, outgoing=="g") and GQ (incoming=="g",
  // outgoing=="q") don't have this problem: only one side ever sums over
  // flavor there, multiplied by a flavor-independent single value on the
  // other side (the gluon PDF for QG, the gluon FF for GQ) -- see
  // point_tables.cpp's QG/GQ branches of integrand_xi, both of which key
  // off pdf.xf() only for the quark side. Summing first and multiplying
  // after is then equivalent to summing the per-flavor products, so no
  // cross-flavor mixing occurs and these two channels are fine as-is.
  if(is_multiflavor_incoming(rp.incoming) && rp.outgoing=="q"){
    std::cerr << "Error: --level hadron with a multi-flavor --incoming (\"" << rp.incoming
               << "\") and --outgoing q is not supported -- it would multiply the "
                  "flavor-summed incoming PDF by the flavor-summed outgoing FF, "
                  "producing unphysical cross-flavor terms (e.g. an incoming "
                  "u-quark fragmenting via the d-quark FF) instead of the "
                  "flavor-diagonal sum the QQ channel actually requires. Run once "
                  "per explicit matching flavor (--incoming u --outgoing u, "
                  "--incoming d --outgoing d, --incoming s --outgoing s, and the "
                  "antiquark pairs) and sum the resulting LO/NLO values instead -- "
                  "see README.md, \"Hadron-level flavor sums\"." << std::endl;
    exit(1);
  }

  gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(n_zpoints);

  // Each Gauss-Legendre node below is independent by construction, so the
  // loop is parallelized over i. LHAPDF::PDF (wrapped by PdfSet/FfSet) is
  // not safe to share across threads -- it lazily builds mutable
  // interpolator/extrapolator caches on first evaluation -- so each thread
  // gets its own PdfSet/FfSet, built serially here from the same
  // underlying set. dipole is safe to share (Spline2D::eval uses a
  // thread_local accelerator, see spline_wrappers.cpp).
  const int nthreads = omp_get_max_threads();
  std::vector<std::unique_ptr<PdfSet>> pdf_tls(nthreads);
  std::vector<std::unique_ptr<FfSet>> ff_tls(nthreads);
  for(int t=0; t<nthreads; t++){
    pdf_tls[t] = std::make_unique<PdfSet>(pdf.name());
    ff_tls[t] = std::make_unique<FfSet>(ff.name());
  }

  // Written one slot per i, then summed serially below in index order --
  // avoids both a shared accumulator and any parallel-reduction reordering,
  // so the result matches the original serial summation exactly.
  std::vector<double> LO_i(n_zpoints, 0.0), NLO_i(n_zpoints, 0.0);

  #pragma omp parallel for schedule(dynamic)
  for(int i=0; i<n_zpoints; i++){
    const PdfSet& tpdf = *pdf_tls[omp_get_thread_num()];
    const FfSet& tff = *ff_tls[omp_get_thread_num()];

    double z, w;
    gsl_integration_glfixed_point(zmin, 1, i, &z, &w, table);

    double k = p_h/z;
    double xp = (k/sqrts)*exp(y);    // Eq. 6, at this z's k
    double xg = (k/sqrts)*exp(-y);   // Eq. 5, at this z's k
    double zD = tff.zD(rp, z, rp.mu2_ff); // Note that LHAPDF returns z*D(z,mu^2)
    double D = zD/z; // below we have explicitly 1/z^2
    double jacobian = w/Sq(z);

    // LO only ever needs S(r,Y) (sigma_LO.cpp), never PointTables' NLO
    // coefficient/xi-convolution tables, so it skips PointTables entirely.
    LO_i[i] = jacobian*D*sigma_LO_k(rp, tpdf, dipole, xg, k, xp);

    // NLO does need the full PointTables, built once here and reused for
    // nothing else -- each Gauss-Legendre node is independent by
    // construction, so there's nothing to cache across them.
    PointTables tables(rp, tpdf, dipole, xp, xg, k);
    NLO_i[i] = jacobian*D*sigma_NLO_k(rp, tpdf, tables, k, xp);
  }

  HadronSigma result{0, 0};
  for(int i=0; i<n_zpoints; i++){
    result.LO += LO_i[i];
    result.NLO += NLO_i[i];
  }

  gsl_integration_glfixed_table_free(table);
  return result;
}
