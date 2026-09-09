#include "init_interp.hpp"

#include "Sr_interp_1D.hpp"
#include "params.hpp"
#include "xi_int.hpp"
#include "xi_interp.hpp"

#include <cmath>

using namespace params;

void init_interp(double xp, double xg, double k){
  //if(with_CF || with_Nc || with_xi1 || with_gl || with_gq || with_gg){
    init_xi_interp(xg);
    init_xi_int_interp(xp,xg,k);
  //}
    init_Sr_interp_1D(log(1/xg));
}


void clear_interp(){
  if(with_CF || with_Nc || with_xi1 || with_gl || with_gq || with_gg){
    clear_xi_interp();
    clear_xi_int_interp();
  }
  clear_Sr_interp_1D();
}
