#pragma once

#include "params.hpp"
#include "spline_wrappers.hpp"

#include <string>
#include <vector>

// Initial-condition parameters read from a BK solution file's header
// comment (the "# Initial condition: MV model, Q_s0^2 = ..., \gamma = ...,
// coefficient of E inside Log is ..., x0=..., \Lambda_QCD = ..." and
// "# Nc=..., Nf=..." lines). Shared by read_initial_condition_header()
// below (used by params::make_run_parameters() to seed
// RunParameters::Qs02/gamm/ec from the proton file) and by
// DipoleAmplitude::load_grid (which records the same fields for whichever
// grid file it actually loads, via its Qs02()/gamma()/ec()/Nf() getters).
struct InitialConditionParams{
  double Qs02 = 0, gamma = 0, ec = 0;
  int Nf = 0;
};

// Reads just the initial-condition header out of a BK solution file --
// not the evolution grid itself. Exits(1) if the file can't be opened, or
// if Q_s0^2/\gamma/the E-log coefficient aren't found in its header (rather
// than silently returning zeros, which would make Sr_0() wrong without
// warning: pow(0,0)==1).
InitialConditionParams read_initial_condition_header(const std::string& filename);

// The BK-evolved dipole amplitude S(r,Y) (arXiv:2310.06640 sec. III):
// reads the BK-solution grid selected by rp.col/rp.b
// (rp.bk_proton_file/rp.bk_nucleus_prefix, defaulting to
// params::bksolpp/bksolpA unless overridden on the CLI) and
// bicubic-interpolates it over (r,Y). Built once per run in main(), after
// the command line has been parsed into a RunParameters.
class DipoleAmplitude{
public:
  explicit DipoleAmplitude(const params::RunParameters& rp);
  DipoleAmplitude(const DipoleAmplitude&) = delete;
  DipoleAmplitude& operator=(const DipoleAmplitude&) = delete;

  // S(r,Y) via bicubic interpolation over the full grid.
  double S(double r, double Y) const;

  // ln(r) integration bounds derived from the grid's r-range, used by the
  // NLO coefficient functions (nlo_coefficients.hpp).
  double min_ln_r() const { return minlnr_; }
  double max_ln_r() const { return maxlnr_; }

  double min_r() const { return minr_; }
  double max_r() const { return maxr_; }
  double x0() const { return x0_; }

  // Initial-condition parameters read from the BK solution file's header
  // comment (the "# Initial condition: MV model, Q_s0^2 = ..., \gamma = ...,
  // coefficient of E inside Log is ..., x0=..., \Lambda_QCD = ..." and
  // "# Nc=..., Nf=..." lines) -- these record what the solver was actually
  // run with, for whichever grid file this instance loaded (rp.col=="pA"
  // may differ from the proton file RunParameters::Qs02/gamm/ec are read
  // from -- see params.hpp).
  double Qs02() const { return qs02_; }
  double gamma() const { return gamma_; }
  double ec() const { return ec_; }
  int Nf() const { return nf_; }

  // The grid's r/Y sampling -- point_tables.cpp builds its own per-r
  // tables (the NLO coefficient functions and their xi-convolution) on
  // this same r sampling.
  int r_points() const { return rpoints_; }
  const std::vector<double>& r_values() const { return rvals_; }
  int y_points() const { return ypoints_; }
  const std::vector<double>& y_values() const { return yvals_; }

private:
  double minr_ = 0, maxr_ = 0, r_mult_ = 0, x0_ = 0;
  double minlnr_ = 0, maxlnr_ = 0;
  int rpoints_ = 0, ypoints_ = 0;
  std::vector<double> rvals_, yvals_;

  // See the Qs02()/gamma()/ec()/Nf() getters above.
  double qs02_ = 0, gamma_ = 0, ec_ = 0;
  int nf_ = 0;

  Spline2D spline_;

  // Reads the BK solution file into Sr[iY][ir] = S(Y,r) and yvals, and
  // sets minr_/r_mult_/rpoints_/x0_ from its header.
  void load_grid(const params::RunParameters& rp,
                  std::vector<std::vector<double>>& Sr,
                  std::vector<double>& yvals_tmp);
};

// Analytic (McLerran-Venugopalan/GBW-like) initial-condition dipole
// amplitude S(r) at the BK evolution starting rapidity X0: the
// dipole-proton amplitude of Eq. (13), or its dipole-nucleus generalization
// via the optical Glauber model, Eq. (14) (rp.col=="pA"; rp.TA is T_A(b)
// from Eq. 14). Used directly by the NLO cross section (sigma_NLO.cpp)
// rather than through the interpolated DipoleAmplitude grid -- a pure
// function of rp, so it needs no DipoleAmplitude instance.
double Sr_0(const params::RunParameters& rp, double r);
