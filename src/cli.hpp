#pragma once

#include <string>

// Command-line interface for the executable (see main.cpp). Two forms are
// accepted:
//
//  - Named-flag form (preferred): every parameter is given as --name value,
//    e.g. `--incoming q --outgoing g`, in any order:
//      --zmin --zmax --zstep --col --b --incoming --outgoing --rc --p --muratio
//    plus optional flags that override params.hpp compile-time constants
//    for this run only:
//      --sqrts (default params::SQRTS), --y (default params::yh)
//      --bk-proton (default params::bksolpp): BK solution file used for
//        --col pp
//      --bk-nucleus (default params::bksolpA): BK solution filename
//        *prefix* used for --col pA -- the code appends the digits of --b
//        directly (no separator is inserted), so a prefix meant to read as
//        "..._<b>" must already end in "_", e.g. --bk-nucleus ./dir/Pb_b_
//      --sigma02 <value> [mb] (default params::sigma0/2): sigma0/2, i.e.
//        half of Sr_0()'s pA-branch normalization (RunParameters::sigma0
//        is stored as 2x this value, matching params::sigma0). The value
//        is taken to already be in GeV^-2 unless a literal trailing "mb"
//        token follows it, in which case it's first converted from
//        millibarns (e.g. --sigma02 47.2 mb).
//      --level <parton|hadron> (default "parton"): parton keeps the
//        existing behavior (main.cpp's z-scan of the bare parton-level
//        cross section, sec. IV). hadron instead convolves that
//        parton-level result with a fragmentation function (sec. V) --
//        --p is then the *hadron's* transverse momentum p_h, --zmin is the
//        lower bound of the (internal, continuous) z=p_h/k integration
//        (its upper bound is always 1), and --zmax/--zstep are accepted
//        but unused. See ff_set.hpp/sigma_hadron.hpp.
//      --ff-set <name> (default "NNFF10_PIsum_nlo"): the LHAPDF
//        fragmentation-function set used when --level hadron; ignored
//        for --level parton.
//      --z-points <n> (default 16, ignored for --level parton): number of
//        nodes in the fixed-order Gauss-Legendre rule used for --level
//        hadron's z integral (sigma_hadron.hpp) -- deterministic and
//        exactly n NLO-cost evaluations, unlike an adaptive quadrature, so
//        this directly controls that mode's runtime (roughly n times the
//        cost of one parton-level NLO point). 16 was picked as a ~1%-or-
//        better default; see the validation note in sigma_hadron.cpp.
//
//  - Legacy positional form, kept for backward compatibility with existing
//    scripts (Script.sh, Script_oberon.sh, tests/regression/cases/*/args):
//      zmin zmax zstep col b incoming outgoing rc p muratio
//    --sqrts/--y have no positional equivalent; params::SQRTS/params::yh are
//    used as-is.
//
// Which form is used is decided by the first token: if it starts with
// "--", the named form is parsed; otherwise the legacy positional form is
// assumed.
namespace cli{

struct Args{
  double zmin, zmax, zstep;
  std::string col;
  double b;
  std::string incoming, outgoing;
  std::string rc;
  double p;
  double muratio;
  double sqrts;   // params::SQRTS unless overridden by --sqrts
  double y;       // params::yh unless overridden by --y
  std::string bk_proton;   // empty unless --bk-proton given (params::make_run_parameters then keeps params::bksolpp)
  std::string bk_nucleus;  // empty unless --bk-nucleus given (params::make_run_parameters then keeps params::bksolpA)
  double sigma0;  // params::sigma0 unless overridden by --sigma02 (already doubled from that flag's sigma0/2 input)
  std::string level;   // "parton" (default) or "hadron"
  std::string ff_set;  // fragmentation-function set name, used only for level=="hadron"
  int z_points;        // Gauss-Legendre node count for level=="hadron"'s z integral (default 16)
};

// Parses argc/argv as passed to main() (argv[0] is the program name).
// Prints a usage message to stderr and exits(1) on any error: an unknown
// flag, a missing required flag/value, the wrong number of positional
// arguments, or an unparsable number.
Args parse(int argc, char* argv[]);

}
