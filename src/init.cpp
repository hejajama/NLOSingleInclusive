#include "init.hpp"

#include "Sr_interp_2D.hpp"
#include "bksol_nlodisfit.hpp"
#include "params.hpp"
#include "utils.hpp"

#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_errno.h>

using namespace params;

void init(const RunParameters& rp){
  LHAPDF::setVerbosity(0);
  pdf::lhpdf=LHAPDF::mkPDF(pdfname,0);

  init_bksol(rp);
  init_Sr_interp_2D();

  gsl_set_error_handler(&gsl_error_handler);
}
