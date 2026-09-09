#pragma once

#include <string>

// Choice of running-coupling prescription for alpha_s(r) in the NLO impact
// factor (arXiv:2310.06640). Selected on the command line (see main.cpp).
enum running_types {FIXED, MOM, PARENT, DAUGHTER, SMALLEST, MIXED, MIXEDBD};

// Run configuration: physical constants, the choice of dipole-amplitude fit
// to use, kinematics, and the mutable state derived from the command-line
// arguments (col, b, incoming, outgoing, alpha_s_running, mu2, TA, ...).
//
// This is still the original flat, partly-mutable global config carried
// over unchanged from the pre-refactor main.cpp — only split into
// declaration (here) + definition (params.cpp) so every translation unit
// can see it without re-defining it. Turning it into a proper immutable
// RunParameters value type is a later refactor stage.
namespace params{
  // First Pb BK solution
  //const string bksolpA = "./Pb/Pb_smallest_b_";
  //const string bksolpp = "./Pb/proton_resumbk_fit1_smallest";

  // KCBK fit1 solution
  //const string bksolpA = "./kcbk_fit1_sigma0_18.81mb/Pb_smallest_b_";
  //const string bksolpp = "./kcbk_fit1_sigma0_18.81mb/p_smallest";

  // KCBK first  bal+sd solution
  //const string bksolpp = "./dipole-kcbk-hera-bal+sd-4.61.dip";

  // KCBK fit3 solution
  extern const std::string bksolpp;
  extern const std::string bksolpA; //"./KCBK_fit_3/Pb_b_";

  // KCBK bal+sd solution
  //const string bksolpp = "./KCBK_fit_5/proton.dat";
  //const string bksolpA = "./KCBK_fit_5/Pb_b_";

  // KCBK parent solution
  //const string bksolpp = "./KCBK_fit_1/proton.dat";
  //const string bksolpA = "./KCBK_fit_1/Pb_b_";

  // ResumBK bal+sd solution
  //const string bksolpp = "./ResumBK_fit_5/proton.dat";
  //const string bksolpA = "./ResumBK_fit_5/Pb_b_";

  // ResumBK parent solution
  //const string bksolpp = "./ResumBK_fit_1/proton";
  //const string bksolpA = "./ResumBK_fit_1/Pb_b_";

  // ResumBK first parent solution
  //const string bksolpp = "./dipole-resumbk-hera-parent-4.61.dip";

  // TBK parent solution
  //const string bksolpp = "./TBK_fit_1/proton.dat";
  //const string bksolpA = "./TBK_fit_1/Pb_b_";

  //const string bksol="/Users/tawabyx/Documents/GitHub/nlodisfit/data/dipole_resumbk_sdrc.dip";
  //const string bksol = "/Users/tawabyx/Documents/code_NLO_sinc/Au/pA_impactb0.000000";
  //const string bksol="bksol/Sgrid_31.dat";
  //AmplitudeLib N(bksol);

  // initial condition
  //const double Qs02=0.2;
  //const double Qs02=0.0964;  //First Pb run and ResumBK parent solution
  //const double Qs02=0.0833;      // KCBK fit1 and parent solution
  extern const double Qs02;      // KCBK fit3
  //const double Qs02=0.0905;	 // KCBK bal+sd solution
  //const double Qs02=0.0950;      // ResumBK bal+sd solution
  //const double Qs02=0.0917;     // TBK parent (fit1) solution
  extern const double ec;
  // running coupling
  //const running_types alpha_s_running=FIXED;
  extern const double alpha_s_fixed; // (alpha_bar=0.2)
  extern const double alpha_s_freeze;
  // choice of terms
  extern bool with_Nc;
  extern bool with_CF;
  extern bool with_gl;
  extern bool with_gq;
  extern bool with_gg;
  /*const bool with_Nc = true;
  const bool with_CF = true;*/
  extern const bool with_xi1;  // True iff one wants to use the subtracted scheme (still not exactly CXY, see the paper)
  /*const bool with_gl = false;    // Include q->g terms
  const bool with_gq = false;   // Include g->q terms
  const bool with_gg = false;   // Include g->g terms */
  // kinematics
  //const double SQRTS=500;    // Forward RHIC (GeV)
  //const double SQRTS=5020;     // LHC (GeV)
  //const double SQRTS=5000;     // LHCb pp (GeV)
  extern const double SQRTS;     // New LHCb (GeV)
  extern const double yh;
  // PDFs
  extern const std::string pdfname;
  extern double mu2;
  // constants
  extern double minlnr, maxlnr;
  extern const int Nc, Nf;

  // const double CF=4./3.;
  extern const double CF;     // To be consistent with the large-Nc limit taken elsewhere

  extern const double beta0;
  extern const double c0;
  extern const double LambdaQCD;
  extern const double alpha_s_mu_0;

  //const int Anucleus = 197;   // Henri's Au data
  extern const int Anucleus;   // Heikki's Pb data

  extern const double RA;
  extern const double WSd;
  extern const double sigma_inel;
  extern const double gamm;
  //const double sigma0 = 50.2628683108;   // Henri's Au data
  //const double sigma0 = 96.6153078;        // ResumBK and old KCBK data
  extern const double sigma0;          // KCBK fit3
  //const double sigma0 = 44.5837784;      // KCBK bal+sd
  //const double sigma0 = 50.0283412;        // KCBK parent
  //const double sigma0 = 39.3446708;        // ResumBK parent
  //const double sigma0 = 40.320583;      // ResumBK bal+sd
  //const double sigma0 = 31.7941922;        // TBK parent (fit1)
  extern double TA;

  extern const std::string TAfile;     // Pb

  // intde parameters
  extern const double epsrel_intde;
  // gsl parameters
  extern const int gsl_maxpoints;
  extern const double epsabs_gsl, epsrel_gsl;

  // input parameters
  extern std::string col;
  extern double b;
  extern double p;
  extern std::string incoming;
  extern std::string outgoing;
  extern running_types alpha_s_running;
}
