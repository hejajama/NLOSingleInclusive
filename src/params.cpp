#include "params.hpp"

#include "gsl/gsl_math.h"
#include <cmath>

namespace params{
  const std::string bksolpp = "./KCBK_fit_3/proton.dat";
  const std::string bksolpA = "./KCBK_fit_3/proton.dat";

  const double Qs02=0.0680;
  const double ec=1.;
  const double alpha_s_fixed=0.2*M_PI/3.;
  const double alpha_s_freeze=0.7;

  bool with_Nc = false;
  bool with_CF = false;
  bool with_gl = false;
  bool with_gq = false;
  bool with_gg = false;
  const bool with_xi1 = false;

  const double SQRTS=8160;
  const double yh=3;

  const std::string pdfname="MSTW2008nlo90cl";
  double mu2;

  double minlnr, maxlnr;
  const int Nc=3, Nf=3;

  const double CF=3./2.;

  const double beta0=(11.*Nc-2.*Nf)/3.;
  const double c0=2*exp(-M_EULER);
  const double LambdaQCD=0.241;
  const double alpha_s_mu_0=exp((2*M_PI)/(beta0*alpha_s_freeze));

  const int Anucleus = 208;

  const double RA = (1.12 * pow(Anucleus, 1/3)) + (0.86 * pow(Anucleus, -1/3));
  const double WSd = 0.54;
  const double sigma_inel = 179.7733;
  const double gamm = 1.21;
  const double sigma0 = 94.4580282;
  double TA;

  const std::string TAfile = "./TAvalues_Pb_Heikki.dat";

  const double epsrel_intde=1e-4;
  const int gsl_maxpoints=1000;
  const double epsabs_gsl=1e-8, epsrel_gsl=1e-6;

  std::string col;
  double b;
  double p;
  std::string incoming;
  std::string outgoing;
  running_types alpha_s_running;
}
