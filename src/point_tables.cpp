#include "point_tables.hpp"

#include "common.hpp"
#include "nlo_coefficients.hpp"
#include "pdf_set.hpp"
#include "running_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

using namespace std;
using namespace params;

namespace {
  // I1 (Eq. 10c) and H1 (Eq. 12a): unlike the rest of the coefficient
  // functions (nlo_coefficients.hpp), these need no transverse integral --
  // they're plain algebra in S(r,Y), so they don't need NLOCoefficients'
  // GSL workspaces and live here as free functions instead. Same symbols
  // as the paper. Local to the xi-convolution -- nothing outside
  // build_xi_convolution()/integrand_xi needs them.
  double I1(const RunParameters& rp, const DipoleAmplitude& dipole, double r, double y){
    double res=(dipole.S(r,y)*(2*log(c0/r)-log(rp.mu2)))/(2*M_PI);
    if(is_position_space_alpha_s(rp.alpha_s_running)){
      return alpha_s_pos(r)*res;
    }else{
      return res;
    }
  }

  double H1(const RunParameters& rp, const DipoleAmplitude& dipole, double r, double y){
    double res=(Sq(dipole.S(r,y))*(2*log(c0/r)-log(rp.mu2)))/M_PI;
    if(is_position_space_alpha_s(rp.alpha_s_running)){
      return alpha_s_pos(r)*res;
    }else{
      return res;
    }
  }

  // Context threaded through integrand_xi's GSL callback.
  struct IntegrandXiContext{
    const RunParameters *rp;
    const PdfSet *pdf;
    const DipoleAmplitude *dipole;
    const PointTables *tables;
    double r, xp, xg, k;
  };
}


PointTables::PointTables(const RunParameters& rp, const PdfSet& pdf,
                          const DipoleAmplitude& dipole, double xp, double xg, double k)
  : minr_(dipole.min_r()), maxr_(dipole.max_r()),
    y_floor_(std::log(1/dipole.x0())),
    outer_slice_(dipole, std::log(1/xg))
{
  build_coefficient_tables(rp, dipole, xg);
  build_xi_convolution(rp, pdf, dipole, xp, xg, k);
}


double PointTables::I2(double r, double y) const{ return coeff_splines_[kI2].eval(r, max(y, y_floor_)); }
double PointTables::J(double r, double y) const{ return coeff_splines_[kJ].eval(r, max(y, y_floor_)); }
double PointTables::K1(double r, double y) const{ return coeff_splines_[kK1].eval(r, max(y, y_floor_)); }
double PointTables::H2(double r, double y) const{ return coeff_splines_[kH2].eval(r, max(y, y_floor_)); }
double PointTables::H3(double r, double y) const{ return coeff_splines_[kH3].eval(r, max(y, y_floor_)); }
double PointTables::H4(double r, double y) const{ return coeff_splines_[kH4].eval(r, max(y, y_floor_)); }
double PointTables::K2(double r, double y) const{ return coeff_splines_[kK2].eval(r, max(y, y_floor_)); }
double PointTables::Jv(double r, double y) const{ return coeff_splines_[kJv].eval(r, max(y, y_floor_)); }
double PointTables::Jv2(double r, double y) const{ return coeff_splines_[kJv2].eval(r, max(y, y_floor_)); }
double PointTables::JJv_xi1(double r, double y) const{ return coeff_splines_[kJJv_xi1].eval(r, max(y, y_floor_)); }


double PointTables::xi_convolution(double r) const{
  if(r<minr_) r=minr_;
  if(r>maxr_) return 0;
  return xi_conv_spline_.eval(r);
}


double PointTables::dipole_slice(double r) const{
  return outer_slice_(r);
}


void PointTables::build_coefficient_tables(const RunParameters& rp, const DipoleAmplitude& dipole, double xg){
  double yg=log(1/xg);

  int ypoints = dipole.y_points();
  const vector<double> &yvals = dipole.y_values();

  int realypoints=0;
  for(int i=0; i<ypoints; i++){
    realypoints++;
    if(1-(yvals[i]+y_floor_)/yg<1e-5) break;
  }

  vector<double> yvals_tmp(realypoints);
  for(int i=0; i<realypoints; i++){
    yvals_tmp[i]=min(yvals[i]+y_floor_,yg);
  }

  int rpoints = dipole.r_points();
  const vector<double> &rvals = dipole.r_values();
  for(auto &spline : coeff_splines_) spline.init(rpoints, realypoints);

  for(int i=0; i<realypoints; i++){
    double y=yvals_tmp[i];
    double xi=1-exp(y-yg);

    DipoleAmplitude1DSlice sr1d(dipole, y);
    NLOCoefficients nlo(dipole);

    for(int j=0; j<rpoints; j++){
      double r=rvals[j];
      double I2_tmp=(rp.channel==Channel::QQ ? nlo.I2(rp,sr1d,r,xi) : 0);
      double J_tmp=(rp.channel==Channel::QQ ? nlo.J(rp,sr1d,r,xi) : 0);
      double K1_tmp=(rp.channel==Channel::QG ? nlo.K1(rp,sr1d,r,xi) : 0);
      double H2_tmp=(rp.channel==Channel::GG ? nlo.H2(rp,sr1d,r,xi) : 0);
      double H3_tmp=(rp.channel==Channel::GG ? nlo.H3(rp,sr1d,r,xi) : 0);
      double H4_tmp=(rp.channel==Channel::GG ? nlo.H4(rp,sr1d,r,xi) : 0);
      double K2_tmp=(rp.channel==Channel::GQ ? nlo.K2(rp,sr1d,r,xi) : 0);
      double Jv_tmp=(rp.channel==Channel::QQ ? nlo.Jv(rp,sr1d,r,xi) : 0);
      double Jv2_tmp=(rp.channel==Channel::QQ ? nlo.Jv2(rp,sr1d,r,xi) : 0);
      double JJv_xi1_tmp=(with_xi1 ? nlo.JJv_xi1(rp,sr1d,r) : 0);

      coeff_splines_[kI2].set(j,i,I2_tmp);
      coeff_splines_[kJ].set(j,i,J_tmp);
      coeff_splines_[kK1].set(j,i,K1_tmp);
      coeff_splines_[kH2].set(j,i,H2_tmp);
      coeff_splines_[kH3].set(j,i,H3_tmp);
      coeff_splines_[kH4].set(j,i,H4_tmp);
      coeff_splines_[kK2].set(j,i,K2_tmp);
      coeff_splines_[kJv].set(j,i,Jv_tmp);
      coeff_splines_[kJv2].set(j,i,Jv2_tmp);
      coeff_splines_[kJJv_xi1].set(j,i,JJv_xi1_tmp);
    }
  }

  for(auto &spline : coeff_splines_) spline.build(rvals.data(), yvals_tmp.data());
}


double PointTables::integrand_xi(double xi, void *userdata){
  const IntegrandXiContext &ctx = *static_cast<IntegrandXiContext*>(userdata);
  const RunParameters &rp = *ctx.rp;
  const PdfSet &pdf = *ctx.pdf;
  const DipoleAmplitude &dipole = *ctx.dipole;
  const PointTables &tables = *ctx.tables;
  double r=ctx.r;
  double xp=ctx.xp;
  double xg=ctx.xg;
  double k=ctx.k;
  //
  double X=xg/(1-xi);
  double y=log(1/X); if(y < 0) y = 0;
  //
  double xi2=Sq(xi);
  double xq_xi=(xi>xp ? pdf.xf(rp,xp/xi,rp.mu2) : 0);
  double xq_xi1=pdf.xf(rp,xp,rp.mu2);
  double xg_xi=(xi>xp ? pdf.xf(rp,xp/xi,rp.mu2) : 0);
  double xg_xi1=pdf.xf(rp,xp,rp.mu2);
  //
  double res=0;
  // qq channel, Eq. (9): the first Channel::QQ block is the I1/I2
  // (C_F-scaling) part and the second is the J/Jv (N_c-scaling) part of the
  // same equation -- the large-Nc decomposition splits Eq. (9) across the
  // two (both execute together, since both are gated on Channel::QQ; kept
  // as two blocks, matching the original with_CF/with_Nc split, rather than
  // merged into one).
  if(rp.channel==Channel::QQ){
    double real_CF=0;
    double virt_CF=-M_1_PI*(2*log(k) - log(rp.mu2) + 2*log(1-xi))*dipole.S(r,y);
    if(is_position_space_alpha_s(rp.alpha_s_running)){
        virt_CF*=alpha_s_pos(r);
        real_CF+=(xi>xp ? I1(rp,dipole,r,y)+I1(rp,dipole,xi*r,y) - 4*tables.I2(r,y) : 0);
    }
    else if(is_mixed_alpha_s(rp.alpha_s_running)){
        virt_CF*=alpha_s_mom(k);
        real_CF+=(xi>xp ? (I1(rp,dipole,r,y)+I1(rp,dipole,xi*r,y) - 4*tables.I2(r,y))*alpha_s_mom(k) : 0);
    }
    else{
        real_CF+=(xi>xp ? I1(rp,dipole,r,y) + I1(rp,dipole,xi*r,y) - 4*tables.I2(r,y) : 0);
    }
    res+=(CF*(1+xi2)*(real_CF*xq_xi+virt_CF*xq_xi1))/(1-xi);
  }
  if(rp.channel==Channel::QQ){
    // qq channel, Eq. (9) -- N_c-scaling part; see the Channel::QQ block above.
    double real_Nc=(xi>xp ? tables.J(r,y) : 0);
    double virt_Nc=-tables.Jv(r,y);
    double resNc=real_Nc*xq_xi+virt_Nc*xq_xi1;
    res+=(Nc*(1+xi2)*resNc)/(1-xi);
  }
  // qg channel, Eq. (11b): [(1/2)I1(xi*r) + (1/4)H1(r) - (1/4)K1(r)].
  if(rp.channel==Channel::QG){
    double resgl1 = (xi>xp ? 0.25*H1(rp,dipole,r,y) + 0.5*I1(rp,dipole,xi*r,y) : 0);
    double resgl2 = (xi>xp ? -0.25*tables.K1(r,y) : 0);
    if(rp.alpha_s_running==PARENT){
        resgl2*=alpha_s_pos(r);
    }
    else if(is_mixed_alpha_s(rp.alpha_s_running)){
        resgl1*=alpha_s_pos(k);
        resgl2*=alpha_s_pos(k);
    }
    res += ((resgl1+resgl2)*Nc*xq_xi*(1+Sq(1-xi)))/xi;
}
  // gq channel, Eq. (11c): [(1/2)I1(r) + (1/4)H1(xi*r) - (1/4)K2(r)].
  if(rp.channel==Channel::GQ){
    double resgq1 = (xi>xp ? 0.25*H1(rp,dipole,xi*r,y) + 0.5*I1(rp,dipole,r,y) : 0);
    double resgq2 = (xi>xp ? -0.25*tables.K2(r,y) : 0);
    if(is_mixed_alpha_s(rp.alpha_s_running)){
        resgq1 *= alpha_s_mom(k);
        if(rp.alpha_s_running == MIXEDBD){
            resgq2 *= alpha_s_mom(k);
        }
    }
    double resgq = resgq1 + resgq2;
    //resgq *= Nf;               // Multiply by Nf to normalize. Technically, no multiplication but sum over flavors with fragmentation fn.
    res += (resgq*xg_xi*(Sq(xi)+Sq(1-xi)));
  }
  // gg channel, Eq. (11a). gg2/gg3 are the two square-bracket terms after
  // the [H1+H1-H2] one; gg3 (Nf-weighted below) is built from tables.H4
  // (Eq. 12d).
  if(rp.channel==Channel::GG){
    double gg11 = (xi>xp ? H1(rp,dipole,r,y) + H1(rp,dipole,xi*r,y) : 0);
    double gg12 = (xi>xp ? - tables.H2(r,y) : 0);
    double gg2 = tables.H3(r,y) - M_1_PI*(2*log(k)-log(rp.mu2)+2*log(1-xi))*Sq(dipole.S(r,y));
    double gg3 = tables.H4(r,y) - 0.5*M_1_PI*(2*log(k)-log(rp.mu2)+2*log(1-xi))*Sq(dipole.S(r,y));
    if(rp.alpha_s_running == PARENT){
        gg2 *= alpha_s_pos(r);
        gg3 *= alpha_s_pos(r);
    }
    else if(is_mixed_alpha_s(rp.alpha_s_running)){
        gg2 *= alpha_s_mom(k);
        gg3 *= alpha_s_mom(k);
        gg11 *= alpha_s_mom(k);
        if (rp.alpha_s_running == MIXEDBD){
            gg12 *= alpha_s_mom(k);
        }
    }
    double resgg = (Nc*xg_xi*(gg11+gg12)*Sq(1-(xi*(1-xi)))/(xi*(1-xi))) + (Nc*xg_xi1*gg2*((2*xi/(1-xi))+(xi*(1-xi))))  + (Nf*xg_xi1*gg3*(Sq(xi)+Sq(1-xi)));
    res += resgg;
  }
  // Never taken: params::with_xi1 is a compile-time false. See
  // nlo_coefficients.hpp for the (unverified) provenance of this term.
  if(with_xi1){
    double sub=Nc*2*tables.JJv_xi1(r,y)*xq_xi1 / (1-xi);
    if(rp.alpha_s_running==PARENT) sub*=alpha_s_pos(r);
    res+=(rp.channel==Channel::QQ ? -sub : sub);
  }
  //res/=1-xi;
  if(gsl_finite(res)==1){
    return res;
  }else{
    return 0;
  }
}


void PointTables::build_xi_convolution(const RunParameters& rp, const PdfSet& pdf,
                                        const DipoleAmplitude& dipole, double xp, double xg, double k){
  int rpoints = dipole.r_points();
  const vector<double> &rvals = dipole.r_values();
  vector<double> intvals(rpoints);

  for(int i=0; i<rpoints; i++){
    double r=rvals[i];
    IntegrandXiContext ctx{&rp, &pdf, &dipole, this, r, xp, xg, k};
    gsl_function F;
    F.function=&PointTables::integrand_xi;
    F.params=&ctx;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(gsl_maxpoints);
    double result, error;
    gsl_integration_qags(&F,0,1-xg,0,epsrel_gsl,gsl_maxpoints,w,&result,&error);
    gsl_integration_workspace_free(w);
    intvals[i]=result;
  }

  xi_conv_spline_.build(rvals.data(), intvals.data(), rpoints);
}
