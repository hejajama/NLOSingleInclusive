#include "init_interp.hpp"

#include "Sr_interp_1D.hpp"
#include "params.hpp"
#include "xi_int.hpp"
#include "xi_interp.hpp"

#include <cmath>

using namespace params;

void init_interp(const RunParameters& rp, double xp, double xg, double k){
  //if(rp.with_CF || rp.with_Nc || with_xi1 || rp.with_gl || rp.with_gq || rp.with_gg){
    init_xi_interp(rp, xg);
    init_xi_int_interp(rp, xp,xg,k);
  //}
    init_Sr_interp_1D(log(1/xg));
}


void clear_interp(const RunParameters& rp){
  if(rp.with_CF || rp.with_Nc || with_xi1 || rp.with_gl || rp.with_gq || rp.with_gg){
    clear_xi_interp();
    clear_xi_int_interp();
  }
  clear_Sr_interp_1D();
}
