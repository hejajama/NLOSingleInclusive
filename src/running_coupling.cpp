#include "running_coupling.hpp"

#include "common.hpp"
#include "params.hpp"

#include <cmath>
#include <gsl/gsl_math.h>

using namespace params;

// Coupling constant runs differently in Henri's dipole. Switch to this one for position space rc.
double alpha_s_pos(double r){
  return (4*M_PI)/(0.2*beta0*log(pow(2.5, 10) + pow(4*1.21/Sq(LambdaQCD*r), 5)));
}


double alpha_s_mom(double Q2){
  return (4*M_PI)/(beta0*log(Q2/Sq(LambdaQCD)));
}


bool is_position_space_alpha_s(running_types rc){
  return rc==PARENT || rc==DAUGHTER || rc==SMALLEST;
}


bool is_mixed_alpha_s(running_types rc){
  return rc==MIXED || rc==MIXEDBD;
}
