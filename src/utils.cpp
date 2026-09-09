#include "utils.hpp"

#include "common.hpp"
#include "params.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

using namespace std;
using namespace params;

PdfSet::PdfSet(const string& name){
  LHAPDF::setVerbosity(0);
  pdf_.reset(LHAPDF::mkPDF(name, 0));
}


double PdfSet::xf(const RunParameters& rp, double x, double Q2) const{
    if(x<0 || x>1){
        cerr << "Bad x= " << x << endl;
        exit(1);
    }

    // 0: gluon
    // 1: d, -1: dbar
    // 2: u, -2: ubar
    // 3: s, -3: sbar
    int particle_flag;
    const string &incoming = rp.incoming;
    if(incoming.compare("g") == 0) particle_flag = 0;
    else if(incoming.compare("d") == 0) particle_flag = 1;
    else if(incoming.compare("u") == 0) particle_flag = 2;
    else if(incoming.compare("s") == 0) particle_flag = 3;
    else if(incoming.compare("dbar") == 0) particle_flag = -1;
    else if(incoming.compare("ubar") == 0) particle_flag = -2;
    else if(incoming.compare("sbar") == 0) particle_flag = -3;

    return pdf_->xfxQ2(particle_flag,x,Q2);
}


// Coupling constant runs differently in Henri's dipole. Switch to this one for position space rc.
double alpha_s_pos(double r){
  return (4*M_PI)/(0.2*beta0*log(pow(2.5, 10) + pow(4*1.21/Sq(LambdaQCD*r), 5)));
}



double alpha_s_mom(double Q2){
  return (4*M_PI)/(beta0*log(Q2/Sq(LambdaQCD)));
}




double str_to_double(string str){
  stringstream stream(str);
  double tmp;
  stream >> tmp;
  return tmp;
}


int str_to_int(string str){
  stringstream stream(str);
  int tmp;
  stream >> tmp;
  return tmp;
}


void gsl_error_handler(const char * reason,
                       const char * file,
                       int line,
                       int error){
  if(error==GSL_EMAXITER || error==GSL_ETOL || error==GSL_EUNDRFLW || error==GSL_EROUND){
    return;
  }else{
    std::cerr << file << ":" << line << ": Error " << error << ": " << reason << std::endl;
    exit(1);
  }
}
