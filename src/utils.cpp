#include "utils.hpp"

#include <cstdlib>
#include <iostream>
#include <gsl/gsl_errno.h>

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
