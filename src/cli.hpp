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
};

// Parses argc/argv as passed to main() (argv[0] is the program name).
// Prints a usage message to stderr and exits(1) on any error: an unknown
// flag, a missing required flag/value, the wrong number of positional
// arguments, or an unparsable number.
Args parse(int argc, char* argv[]);

}
