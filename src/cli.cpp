#include "cli.hpp"
#include "params.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

namespace cli{

namespace{

[[noreturn]] void usage_error(const string& prog, const string& message){
  cerr << "Error: " << message << "\n\n"
       << "Usage (named flags):\n"
       << "  " << prog << " --zmin <d> --zmax <d> --zstep <d> --col <pp|pA> --b <d>\n"
       << "      --incoming <g|q|u|d|s|ubar|dbar|sbar> --outgoing <g|q>\n"
       << "      --rc <fixed|mom|parent|daughter|smallest|mixed|mixedbd>\n"
       << "      --p <d> --muratio <d> [--sqrts <d>] [--y <d>]\n"
       << "      [--bk-proton <file>] [--bk-nucleus <prefix>]\n"
       << "      [--sigma02 <d> [mb]]\n"
       << "      [--level <parton|hadron>] [--ff-set <name>] [--z-points <n>]\n\n"
       << "Usage (legacy positional, kept for backward compatibility):\n"
       << "  " << prog << " zmin zmax zstep col b incoming outgoing rc p muratio\n";
  exit(1);
}

double parse_double(const string& prog, const string& name, const string& value){
  try{
    return stod(value);
  }
  catch(...){
    usage_error(prog, "invalid numeric value for " + name + ": '" + value + "'");
  }
}

int parse_int(const string& prog, const string& name, const string& value){
  try{
    return stoi(value);
  }
  catch(...){
    usage_error(prog, "invalid integer value for " + name + ": '" + value + "'");
  }
}

} // namespace

Args parse(int argc, char* argv[]){
  const string prog = (argc > 0) ? argv[0] : "nlosingleinclusive";
  const vector<string> tokens(argv+1, argv+argc);

  Args args;
  args.sqrts = params::SQRTS;
  args.y = params::yh;
  args.bk_proton = "";
  args.bk_nucleus = "";
  args.sigma0 = params::sigma0;
  args.level = "parton";
  args.ff_set = "NNFF10_PIsum_nlo";
  args.z_points = 16;

  const bool named = !tokens.empty() && tokens[0].rfind("--", 0) == 0;

  if(!named){
    // Legacy positional form: zmin zmax zstep col b incoming outgoing rc p muratio
    if(tokens.size() != 10){
      usage_error(prog, "expected 10 positional arguments (zmin zmax zstep "
                         "col b incoming outgoing rc p muratio), got " +
                         to_string(tokens.size()));
    }
    args.zmin = parse_double(prog, "zmin", tokens[0]);
    args.zmax = parse_double(prog, "zmax", tokens[1]);
    args.zstep = parse_double(prog, "zstep", tokens[2]);
    args.col = tokens[3];
    args.b = parse_double(prog, "b", tokens[4]);
    args.incoming = tokens[5];
    args.outgoing = tokens[6];
    args.rc = tokens[7];
    args.p = parse_double(prog, "p", tokens[8]);
    args.muratio = parse_double(prog, "muratio", tokens[9]);
    return args;
  }

  // Named-flag form: --name value, in any order. --sigma02 is the one
  // variable-arity flag: it optionally takes a second token "mb" right
  // after its number, marking that number as millibarns (to be converted
  // to GeV^-2) rather than already being in GeV^-2.
  map<string,string> flags;
  bool sigma02_is_mb = false;
  for(size_t i=0; i<tokens.size(); i++){
    const string& tok = tokens[i];
    if(tok.rfind("--", 0) != 0){
      usage_error(prog, "expected a --flag, got '" + tok + "'");
    }
    const string name = tok.substr(2);
    if(i+1 >= tokens.size()){
      usage_error(prog, "flag --" + name + " is missing its value");
    }
    flags[name] = tokens[++i];

    if(name == "sigma02" && i+1 < tokens.size()){
      string unit = tokens[i+1];
      transform(unit.begin(), unit.end(), unit.begin(),
                [](unsigned char c){ return std::tolower(c); });
      if(unit == "mb"){
        sigma02_is_mb = true;
        i++;  // consume the unit token too
      }
    }
  }

  auto require = [&](const string& name) -> const string&{
    auto it = flags.find(name);
    if(it == flags.end()){
      usage_error(prog, "missing required flag --" + name);
    }
    return it->second;
  };

  args.zmin = parse_double(prog, "zmin", require("zmin"));
  args.zmax = parse_double(prog, "zmax", require("zmax"));
  args.zstep = parse_double(prog, "zstep", require("zstep"));
  args.col = require("col");
  args.b = parse_double(prog, "b", require("b"));
  args.incoming = require("incoming");
  args.outgoing = require("outgoing");
  args.rc = require("rc");
  args.p = parse_double(prog, "p", require("p"));
  args.muratio = parse_double(prog, "muratio", require("muratio"));
  if(flags.count("sqrts")) args.sqrts = parse_double(prog, "sqrts", flags.at("sqrts"));
  if(flags.count("y")) args.y = parse_double(prog, "y", flags.at("y"));
  if(flags.count("bk-proton")) args.bk_proton = flags.at("bk-proton");
  if(flags.count("bk-nucleus")) args.bk_nucleus = flags.at("bk-nucleus");
  if(flags.count("sigma02")){
    double sigma02 = parse_double(prog, "sigma02", flags.at("sigma02"));
    if(sigma02_is_mb) sigma02 /= params::GeV2_to_mb;  // mb -> GeV^-2
    args.sigma0 = 2.0*sigma02;
  }
  if(flags.count("level")){
    args.level = flags.at("level");
    if(args.level != "parton" && args.level != "hadron"){
      usage_error(prog, "--level must be 'parton' or 'hadron', got '" + args.level + "'");
    }
  }
  if(flags.count("ff-set")) args.ff_set = flags.at("ff-set");
  if(flags.count("z-points")){
    args.z_points = parse_int(prog, "z-points", flags.at("z-points"));
    if(args.z_points < 1){
      usage_error(prog, "--z-points must be at least 1, got " + to_string(args.z_points));
    }
  }

  // Reject unknown flags (typos) rather than silently ignoring them.
  static const vector<string> known = {"zmin", "zmax", "zstep", "col", "b",
    "incoming", "outgoing", "rc", "p", "muratio", "sqrts", "y",
    "bk-proton", "bk-nucleus", "sigma02", "level", "ff-set", "z-points"};
  for(const auto& [name, value] : flags){
    (void)value;
    if(find(known.begin(), known.end(), name) == known.end()){
      usage_error(prog, "unknown flag --" + name);
    }
  }

  return args;
}

}
