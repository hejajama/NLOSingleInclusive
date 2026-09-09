#pragma once

#include <cmath>
#include <string>
#include <gsl/gsl_math.h>

// Choice of running-coupling prescription for alpha_s(r) in the NLO impact
// factor (arXiv:2310.06640). Selected on the command line (see main.cpp).
enum running_types {FIXED, MOM, PARENT, DAUGHTER, SMALLEST, MIXED, MIXEDBD};

// Parton-level channel: which combination of incoming/outgoing partons a
// run computes, derived from the "incoming"/"outgoing" CLI tokens in
// params::make_run_parameters(). Replaces the five with_Nc/with_CF/with_gl/
// with_gq/with_gg booleans RunParameters used to carry (with_Nc and
// with_CF were always set together, for exactly the QQ channel below --
// see docs/PAPER_MAPPING.md, "Channel -> equation -> RunParameters flag",
// for which paper equation each value selects: QQ -> Eq. (9), GG ->
// Eq. (11a), QG -> Eq. (11b) [old with_gl], GQ -> Eq. (11c) [old with_gq]).
enum class Channel { QQ, QG, GQ, GG };

// Physics constants and analysis choices fixed by this particular study and
// BK-fit choice (see the comments below for the alternatives this has been
// swapped between) -- change by editing here and recompiling. These never
// vary within a run, so unlike RunParameters (below) they are plain
// compile-time-ish constants, not something threaded through the call
// graph: C++17 inline variables let every translation unit share exactly
// one definition without a separate params.cpp.
namespace params{
  // First Pb BK solution
  //inline const std::string bksolpA = "./Pb/Pb_smallest_b_";
  //inline const std::string bksolpp = "./Pb/proton_resumbk_fit1_smallest";

  // KCBK fit1 solution
  //inline const std::string bksolpA = "./kcbk_fit1_sigma0_18.81mb/Pb_smallest_b_";
  //inline const std::string bksolpp = "./kcbk_fit1_sigma0_18.81mb/p_smallest";

  // KCBK first  bal+sd solution
  //inline const std::string bksolpp = "./dipole-kcbk-hera-bal+sd-4.61.dip";

  // KCBK fit3 solution
  inline const std::string bksolpp = "./KCBK_fit_3/proton.dat";
  inline const std::string bksolpA = "./KCBK_fit_3/Pb_b_";
  // NOTE: Qs02/gamm/ec used to be hand-set here per BK-fit choice (see git
  // history). They're now read straight out of bksolpp's/bk_proton's own
  // "# Initial condition: ..." header comment -- see
  // RunParameters::Qs02/gamm/ec (populated by make_run_parameters() via
  // read_initial_condition_header(), dipole_amplitude.hpp) -- so switching
  // bksolpp above (or --bk-proton) is enough; nothing to keep in sync here.

  // KCBK bal+sd solution
  //inline const std::string bksolpp = "./KCBK_fit_5/proton.dat";
  //inline const std::string bksolpA = "./KCBK_fit_5/Pb_b_";

  // KCBK parent solution
  //inline const std::string bksolpp = "./KCBK_fit_1/proton.dat";
  //inline const std::string bksolpA = "./KCBK_fit_1/Pb_b_";

  // ResumBK bal+sd solution
  //inline const std::string bksolpp = "./ResumBK_fit_5/proton.dat";
  //inline const std::string bksolpA = "./ResumBK_fit_5/Pb_b_";

  // ResumBK parent solution
  //inline const std::string bksolpp = "./ResumBK_fit_1/proton";
  //inline const std::string bksolpA = "./ResumBK_fit_1/Pb_b_";

  // ResumBK first parent solution
  //inline const std::string bksolpp = "./dipole-resumbk-hera-parent-4.61.dip";

  // TBK parent solution
  //inline const std::string bksolpp = "./TBK_fit_1/proton.dat";
  //inline const std::string bksolpA = "./TBK_fit_1/Pb_b_";

  // running coupling
  //inline constexpr running_types alpha_s_running = FIXED;
  inline constexpr double alpha_s_fixed = 0.2*M_PI/3.; // (alpha_bar=0.2)
  inline constexpr double alpha_s_freeze = 0.7;
  // choice of terms: with_xi1 is a fixed analysis choice; RunParameters::channel
  // is *derived from the command line* (params::RunParameters below), since
  // it selects which parton channel (qq/qg/gq/gg) is being computed.
  inline constexpr bool with_xi1 = false;  // True iff one wants to use the subtracted scheme (still not exactly CXY, see the paper)
  // kinematics
  //inline constexpr double SQRTS = 500;    // Forward RHIC (GeV)
  //inline constexpr double SQRTS = 5020;     // LHC (GeV)
  //inline constexpr double SQRTS = 5000;     // LHCb pp (GeV)
  inline constexpr double SQRTS = 8160;     // New LHCb (GeV)
  inline constexpr double yh = 3;
  // PDFs
  inline const std::string pdfname = "MSTW2008nlo90cl";
  // constants
  inline constexpr int Nc = 3, Nf = 3;

  // inline constexpr double CF = 4./3.;
  inline constexpr double CF = 3./2.;     // To be consistent with the large-Nc limit taken elsewhere

  inline constexpr double beta0 = (11.*Nc-2.*Nf)/3.;
  inline const double c0 = 2*std::exp(-M_EULER); // exp() isn't constexpr, so this can't be `constexpr`
  inline constexpr double LambdaQCD = 0.241;
  inline const double alpha_s_mu_0 = std::exp((2*M_PI)/(beta0*alpha_s_freeze));

  //inline constexpr int Anucleus = 197;   // Henri's Au data
  inline constexpr int Anucleus = 208;   // Heikki's Pb data

  inline const double RA = (1.12 * std::pow(Anucleus, 1/3)) + (0.86 * std::pow(Anucleus, -1/3));
  inline constexpr double WSd = 0.54;
  inline constexpr double sigma_inel = 179.7733;
  //inline constexpr double sigma0 = 50.2628683108;   // Henri's Au data
  //inline constexpr double sigma0 = 96.6153078;        // ResumBK and old KCBK data
  inline constexpr double sigma0 = 94.4580282;          // KCBK fit3
  //inline constexpr double sigma0 = 44.5837784;      // KCBK bal+sd
  //inline constexpr double sigma0 = 50.0283412;        // KCBK parent
  //inline constexpr double sigma0 = 39.3446708;        // ResumBK parent
  //inline constexpr double sigma0 = 40.320583;      // ResumBK bal+sd
  //inline constexpr double sigma0 = 31.7941922;        // TBK parent (fit1)

  inline const std::string TAfile = "./TAvalues_Pb_Heikki.dat";     // Pb

  // intde parameters
  inline constexpr double epsrel_intde = 1e-4;
  // gsl parameters
  inline constexpr int gsl_maxpoints = 1000;
  inline constexpr double epsabs_gsl = 1e-8, epsrel_gsl = 1e-6;

  // Everything above is fixed for the life of the program. Everything below
  // is derived once from the command line (main.cpp) and then read-only for
  // the rest of the run -- built by parse_run_parameters() and passed as
  // `const RunParameters&` through the call graph instead of being reached
  // via a blanket `using namespace params;` the way it used to be.
  struct RunParameters{
    std::string col;
    double b;
    double p;
    std::string incoming;
    std::string outgoing;
    running_types alpha_s_running;
    double mu2;
    double TA;              // params::TAfile row matching b
    Channel channel;
    std::string bk_proton_file;    // BK solution file for col=="pp" (params::bksolpp unless overridden on the CLI)
    std::string bk_nucleus_prefix; // BK solution filename *prefix* for col=="pA" (params::bksolpA unless overridden); load_grid() appends the digits of b directly, so a prefix meant to look like "..._<b>" must already end in "_"

    // Initial-condition parameters for Sr_0() (dipole_amplitude.hpp), read
    // out of bk_proton_file's own header comment by make_run_parameters()
    // (via read_initial_condition_header()) -- always from the *proton*
    // file, even when col=="pA", since Sr_0's nucleus branch is the same
    // proton initial condition generalized through the optical Glauber TA
    // factor, not a separate nucleus Qs02/gamma/ec.
    double Qs02, gamm, ec;
  };

  // Looks up TA for the given impact parameter b in params::TAfile.
  double lookup_TA(double b);

  // Parses the "rc" command-line token ("fixed"/"mom"/"parent"/"daughter"/
  // "mixed"/"mixedbd"/"smallest") into a running_types value. Prints a
  // warning and defaults to FIXED on an unrecognized token, matching the
  // original (silent, zero-initialization-based) fallback exactly.
  running_types parse_alpha_s_running(const std::string& rc);

  // Builds a RunParameters from the already-parsed CLI fields: looks up TA
  // via lookup_TA(b), and derives channel from incoming/outgoing.
  // bk_proton/bk_nucleus override params::bksolpp/params::bksolpA for this
  // run when non-empty (see --bk-proton/--bk-nucleus in cli.hpp); an empty
  // string (the default) keeps the compiled-in params.hpp value.
  RunParameters make_run_parameters(std::string col, double b, double p,
                                     std::string incoming, std::string outgoing,
                                     running_types alpha_s_running, double mu2,
                                     const std::string& bk_proton = "",
                                     const std::string& bk_nucleus = "");
}
