//  Copyright (C) 2018 Yi Pan <ypan1988@gmail.com>
//
//  This program is free software from derived from roptim Package (https://cran.r-project.org/web/packages/roptim/index.html) 
//  under their Licenses
//  https://www.R-project.org/Licenses/


#ifndef FUNCTOR_H_
#define FUNCTOR_H_

#include <RcppArmadillo.h>

namespace roptim {

struct OptStruct {
  bool has_grad_ = false;
  bool has_hess_ = false;
  arma::vec ndeps_;       // tolerances for numerical derivatives
  double fnscale_ = 1.0;  // scaling for objective
  arma::vec parscale_;    // scaling for parameters
  int usebounds_ = 0;
  arma::vec lower_, upper_;
  bool sann_use_custom_function_ = false;
};

class Functor {
public:
  Functor() {  // UpdateOptStruct();
  }
  
  virtual ~Functor() {}
  
  virtual double operator()(const arma::vec &par) = 0;
  virtual void Gradient(const arma::vec &par, arma::vec &grad) {
    ApproximateGradient(par, grad);
  }
  virtual void Hessian(const arma::vec &par, arma::mat &hess) {
    ApproximateHessian(par, hess);
  }
  
  // Returns forward-difference approximation of Gradient
  void ApproximateGradient(const arma::vec &par, arma::vec &grad);
  
  // Returns forward-difference approximation of Hessian
  void ApproximateHessian(const arma::vec &par, arma::mat &hess);
  
  OptStruct os;
};

inline void Functor::ApproximateGradient(const arma::vec &par,
                                         arma::vec &grad) {
  if (os.parscale_.is_empty()) os.parscale_ = arma::ones<arma::vec>(par.size());
  if (os.ndeps_.is_empty())
    os.ndeps_ = arma::ones<arma::vec>(par.size()) * 1e-3;
  
  grad = arma::zeros<arma::vec>(par.size());
  arma::vec x = par % os.parscale_;
  
  if (os.usebounds_ == 0) {
    for (std::size_t i = 0; i != par.size(); ++i) {
      double eps = os.ndeps_(i);
      
      x(i) = (par(i) + eps) * os.parscale_(i);
      double val1 = operator()(x) / os.fnscale_;
      
      x(i) = (par(i) - eps) * os.parscale_(i);
      double val2 = operator()(x) / os.fnscale_;
      
      grad(i) = (val1 - val2) / (2 * eps);
      
      x(i) = par(i) * os.parscale_(i);
    }
  } else {  // use bounds
    for (std::size_t i = 0; i != par.size(); ++i) {
      double epsused = os.ndeps_(i);
      double eps = os.ndeps_(i);
      
      double tmp = par(i) + eps;
      if (tmp > os.upper_(i)) {
        tmp = os.upper_(i);
        epsused = tmp - par(i);
      }
      
      x(i) = tmp * os.parscale_(i);
      double val1 = operator()(x) / os.fnscale_;
      
      tmp = par(i) - eps;
      if (tmp < os.lower_(i)) {
        tmp = os.lower_(i);
        eps = par(i) - tmp;
      }
      
      x(i) = tmp * os.parscale_(i);
      double val2 = operator()(x) / os.fnscale_;
      
      grad(i) = (val1 - val2) / (epsused + eps);
      
      x(i) = par(i) * os.parscale_(i);
    }
  }
}

inline void Functor::ApproximateHessian(const arma::vec &par, arma::mat &hess) {
  if (os.parscale_.is_empty()) os.parscale_ = arma::ones<arma::vec>(par.size());
  if (os.ndeps_.is_empty())
    os.ndeps_ = arma::ones<arma::vec>(par.size()) * 1e-3;
  
  hess = arma::zeros<arma::mat>(par.size(), par.size());
  arma::vec dpar = par / os.parscale_;
  arma::vec df1 = arma::zeros<arma::vec>(par.size());
  arma::vec df2 = arma::zeros<arma::vec>(par.size());
  
  for (std::size_t i = 0; i != par.size(); ++i) {
    double eps = os.ndeps_(i) / os.parscale_(i);
    dpar(i) += eps;
    Gradient(dpar, df1);
    dpar(i) -= 2 * eps;
    Gradient(dpar, df2);
    for (std::size_t j = 0; j != par.size(); ++j)
      hess(i, j) = os.fnscale_ * (df1(j) - df2(j)) /
        (2 * eps * os.parscale_(i) * os.parscale_(j));
    dpar(i) = dpar(i) + eps;
  }
  
  // now symmetrize
  for (std::size_t i = 0; i != par.size(); ++i) {
    for (std::size_t j = 0; j != par.size(); ++j) {
      double tmp = 0.5 * (hess(i, j) + hess(j, i));
      
      hess(i, j) = tmp;
      hess(j, i) = tmp;
    }
  }
}

inline double fminfn(int n, double *x, void *ex) {
  OptStruct os(static_cast<Functor *>(ex)->os);
  
  arma::vec par(x, n);
  par %= os.parscale_;
  return static_cast<Functor *>(ex)->operator()(par) / os.fnscale_;
}

inline void fmingr(int n, double *x, double *gr, void *ex) {
  OptStruct os(static_cast<Functor *>(ex)->os);
  
  arma::vec par(x, n), grad;
  par %= os.parscale_;
  static_cast<Functor *>(ex)->Gradient(par, grad);
  for (auto i = 0; i != n; ++i)
    gr[i] = grad(i) * (os.parscale_(i) / os.fnscale_);
}

}  // namespace roptim

#endif


#ifndef APPLIC_H_
#define APPLIC_H_

// This file is part of R_ext/Applic.h for the function prototype of
// optimization routines. I am not able to include R_ext/Applic.h in roptim.h
// directly since it conflict with RcppArmadillo.h.

#ifdef __cplusplus
extern "C" {
#endif
  
  typedef double optimfn(int, double *, void *);
  typedef void optimgr(int, double *, double *, void *);
  
  void vmmin(int n, double *b, double *Fmin, optimfn fn, optimgr gr, int maxit,
             int trace, int *mask, double abstol, double reltol, int nREPORT,
             void *ex, int *fncount, int *grcount, int *fail);
  void nmmin(int n, double *Bvec, double *X, double *Fmin, optimfn fn, int *fail,
             double abstol, double intol, void *ex, double alpha, double bet,
             double gamm, int trace, int *fncount, int maxit);
  void cgmin(int n, double *Bvec, double *X, double *Fmin, optimfn fn, optimgr gr,
             int *fail, double abstol, double intol, void *ex, int type,
             int trace, int *fncount, int *grcount, int maxit);
  void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
              double *Fmin, optimfn fn, optimgr gr, int *fail, void *ex,
              double factr, double pgtol, int *fncount, int *grcount, int maxit,
              char *msg, int trace, int nREPORT);
  void samin(int n, double *pb, double *yb, optimfn fn, int maxit, int tmax,
             double ti, int trace, void *ex);
  
#ifdef __cplusplus
}
#endif

#endif  // APPLIC_H_









#ifndef SAMIN_H_
#define SAMIN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <R_ext/Print.h>   // Rprintf
#include <R_ext/Random.h>  // random number generation in samin()
#include <Rinternals.h>

#include <RcppArmadillo.h>


namespace roptim {
namespace internal {

static double *vect(int n) { return (double *)R_alloc(n, sizeof(double)); }

static void genptry(int n, double *p, double *ptry, double scale, void *ex) {
  SEXP s, x;
  int i;
  OptStruct OS = static_cast<Functor *>(ex)->os;
  PROTECT_INDEX ipx;
  
  if (OS.sann_use_custom_function_) {
    /* user defined generation of candidate point */
    PROTECT(x = Rf_allocVector(REALSXP, n));
    arma::vec x_copy = arma::zeros<arma::vec>(n);
    for (i = 0; i < n; i++) {
      if (!R_FINITE(p[i])) Rf_error("non-finite value supplied by 'optim'");
      REAL(x)[i] = p[i] * (OS.parscale_(i));
      x_copy(i) = REAL(x)[i];
    }
    arma::vec grad;
    static_cast<Functor *>(ex)->Gradient(x_copy, grad);
    PROTECT_WITH_INDEX(s = Rcpp::wrap(grad), &ipx);
    REPROTECT(s = Rf_coerceVector(s, REALSXP), ipx);
    if (LENGTH(s) != n)
      Rf_error("candidate point in 'optim' evaluated to length %d not %d",
               LENGTH(s), n);
    for (i = 0; i < n; i++) ptry[i] = REAL(s)[i] / (OS.parscale_(i));
    UNPROTECT(2);
  } else { /* default Gaussian Markov kernel */
    for (i = 0; i < n; i++)
      ptry[i] = p[i] + scale * norm_rand(); /* new candidate point */
  }
}

inline void samin(int n, double *pb, double *yb, optimfn fminfn, int maxit,
                  int tmax, double ti, int trace, void *ex)

/* Given a starting point pb[0..n-1], simulated annealing minimization
 is performed on the function fminfn. The starting temperature
 is input as ti. To make sann work silently set trace to zero.
 sann makes in total maxit function evaluations, tmax
 evaluations at each temperature. Returned quantities are pb
 (the location of the minimum), and yb (the minimum value of
 the function func).  Author: Adrian Trapletti
 */
{
  double E1 = 1.7182818; /* exp(1.0)-1.0 */
double big = 1.0e+35;  /*a very large number*/

long j;
int k, its, itdoc;
double t, y, dy, ytry, scale;
double *p, *ptry;

/* Above have: if(trace != 0) trace := REPORT control argument = STEPS */
if (trace < 0) Rf_error("trace, REPORT must be >= 0 (method = \"SANN\")");

if (n == 0) { /* don't even attempt to optimize */
*yb = fminfn(n, pb, ex);
  return;
}
p = vect(n);
ptry = vect(n);
GetRNGstate();
*yb = fminfn(n, pb, ex); /* init best system state pb, *yb */
if (!R_FINITE(*yb)) *yb = big;
for (j = 0; j < n; j++) p[j] = pb[j];
y = *yb; /* init system state p, y */
if (trace) {
  Rprintf("sann objective function values\n");
  Rprintf("initial       value %f\n", *yb);
}
scale = 1.0 / ti;
its = itdoc = 1;
while (its < maxit) {             /* cool down system */
t = ti / log((double)its + E1); /* temperature annealing schedule */
k = 1;
while ((k <= tmax) && (its < maxit)) /* iterate at constant temperature */
{
  genptry(n, p, ptry, scale * t, ex); /* generate new candidate point */
ytry = fminfn(n, ptry, ex);
if (!R_FINITE(ytry)) ytry = big;
dy = ytry - y;
if ((dy <= 0.0) || (unif_rand() < exp(-dy / t))) { /* accept new point? */
for (j = 0; j < n; j++) p[j] = ptry[j];
  y = ytry;     /* update system state p, y */
if (y <= *yb) /* if system state is best, then update best system state
 pb, *yb */
{
  for (j = 0; j < n; j++) pb[j] = p[j];
  *yb = y;
}
}
its++;
k++;
}
if (trace && ((itdoc % trace) == 0))
  Rprintf("iter %8d value %f\n", its - 1, *yb);
itdoc++;
}
if (trace) {
  Rprintf("final         value %f\n", *yb);
  Rprintf("sann stopped after %d iterations\n", its - 1);
}
PutRNGstate();
}

}  // namespace internal
}  // namespace roptim

#endif










#ifndef ROPTIM_H_
#define ROPTIM_H_

#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif




namespace roptim {

template <typename Derived>
class Roptim {
public:
  std::string method_;
  arma::vec lower_, upper_;
  bool hessian_flag_ = false;
  arma::mat hessian_;
  
  arma::vec lower() const { return lower_; }
  arma::vec upper() const { return upper_; }
  
  arma::vec par() const { return par_; }
  double value() const { return val_; }
  int fncount() const { return fncount_; }
  int grcount() const { return grcount_; }
  int convergence() const { return fail_; }
  std::string message() const { return message_; }
  arma::mat hessian() const { return hessian_; }
  
  void print() const;
  
private:
  arma::vec par_;
  double val_ = 0.0;
  int fncount_ = 0;
  int grcount_ = 0;
  int fail_ = 0;
  std::string message_ = "NULL";
  
public:
  struct RoptimControl {
    std::size_t trace = 0;
    double fnscale = 1.0;
    arma::vec parscale;
    arma::vec ndeps;
    std::size_t maxit = 100;
    double abstol = R_NegInf;
    double reltol = sqrt(2.220446e-16);
    double alpha = 1.0;
    double beta = 0.5;
    double gamma = 2.0;
    int REPORT = 10;
    bool warn_1d_NelderMead = true;
    int type = 1;
    int lmm = 5;
    double factr = 1e7;
    double pgtol = 0.0;
    double temp = 10.0;
    int tmax = 10;
  } control;
  
  Roptim(const std::string method = "Nelder-Mead") : method_(method) {
    if (method_ != "Nelder-Mead" && method_ != "BFGS" && method_ != "CG" &&
        method_ != "L-BFGS-B" && method_ != "SANN")
      Rcpp::stop("Roptim::Roptim(): unknown 'method'");
    
    // Sets default value for maxit & REPORT (which depend on method)
    if (method_ == "Nelder-Mead") {
      control.maxit = 500;
    } else if (method_ == "SANN") {
      control.maxit = 10000;
      control.REPORT = 100;
    }
  }
  
  void set_method(const std::string &method) {
    if (method != "Nelder-Mead" && method != "BFGS" && method != "CG" &&
        method != "L-BFGS-B" && method != "SANN")
      Rcpp::stop("Roptim::set_method(): unknown 'method'");
    else
      method_ = method;
    
    // Sets default value for maxit & REPORT (which depend on method)
    if (method_ == "Nelder-Mead") {
      control.maxit = 500;
      control.REPORT = 10;
    } else if (method_ == "SANN") {
      control.maxit = 10000;
      control.REPORT = 100;
    } else {
      control.maxit = 100;
      control.REPORT = 10;
    }
  }
  
  void set_lower(const arma::vec &lower) {
    if (method_ != "L-BFGS-B")
      Rcpp::warning(
        "Roptim::set_lower(): bounds can only be used with method L-BFGS-B");
    method_ = "L-BFGS-B";
    lower_ = lower;
  }
  
  void set_upper(const arma::vec &upper) {
    if (method_ != "L-BFGS-B")
      Rcpp::warning(
        "Roptim::set_upper(): bounds can only be used with method L-BFGS-B");
    method_ = "L-BFGS-B";
    upper_ = upper;
  }
  
  void set_hessian(bool flag) { hessian_flag_ = flag; }
  
  void minimize(Derived &func, arma::vec &par);
};

template <typename Derived>
inline void Roptim<Derived>::print() const {
  par_.t().print(".par()");
  Rcpp::Rcout << "\n.value()\n" << val_ << std::endl;
  Rcpp::Rcout << "\n.fncount()\n" << fncount_ << std::endl;
  
  if (method_ == "Nelder-Mead" || method_ == "SANN")
    Rcpp::Rcout << "\n.grcount()\nNA" << std::endl;
  else
    Rcpp::Rcout << "\n.grcount()\n" << grcount_ << std::endl;
  
  Rcpp::Rcout << "\n.convergence()\n" << fail_ << std::endl;
  Rcpp::Rcout << "\n.message()\n" << message_ << std::endl;
  if (hessian_flag_) hessian_.print("\n.hessian()");
  Rcpp::Rcout << std::endl;
}

template <typename Derived>
inline void Roptim<Derived>::minimize(Derived &func, arma::vec &par) {
  // PART 1: optim()
  
  // Checks if lower and upper bounds is used
  if ((!lower_.is_empty() || !upper_.is_empty()) && method_ != "L-BFGS-B") {
    Rcpp::warning("bounds can only be used with method L-BFGS-B");
    method_ = "L-BFGS-B";
  }
  
  // Sets the parameter size
  std::size_t npar = par.size();
  
  // Sets default value for parscale & ndeps (which depend on npar)
  if (control.parscale.is_empty())
    control.parscale = arma::ones<arma::vec>(npar);
  if (control.ndeps.is_empty())
    control.ndeps = arma::ones<arma::vec>(npar) * 1e-3;
  
  // Checks control variable trace
  if (method_ == "SANN" && control.trace && control.REPORT == 0)
    Rcpp::stop("'trace != 0' needs 'REPORT >= 1'");
  
  // Note that "method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol'
  // and 'abstol'". There is no simple way to detect whether users set new
  // values for 'reltol' and 'abstol'.
  
  // Gives warning of 1-dim optimization by Nelder-Mead
  if (npar == 1 && method_ == "Nelder-Mead" && control.warn_1d_NelderMead)
    Rcpp::warning("one-dimensional optimization by Nelder-Mead is unreliable");
  
  // Sets default value for lower_
  if (method_ == "L-BFGS-B" && lower_.is_empty()) {
    lower_ = arma::zeros<arma::vec>(npar);
    lower_.for_each([](arma::mat::elem_type &val) { val = R_NegInf; });
  }
  // Sets default value for upper_
  if (method_ == "L-BFGS-B" && upper_.is_empty()) {
    upper_ = arma::zeros<arma::vec>(npar);
    upper_.for_each([](arma::mat::elem_type &val) { val = R_PosInf; });
  }
  
  // PART 2: C_optim()
  
  func.os.usebounds_ = 0;
  func.os.fnscale_ = control.fnscale;
  func.os.parscale_ = control.parscale;
  
  if (control.ndeps.size() != npar)
    Rcpp::stop("'ndeps' is of the wrong length");
  else
    func.os.ndeps_ = control.ndeps;
  
  arma::vec dpar = arma::zeros<arma::vec>(npar);
  arma::vec opar = arma::zeros<arma::vec>(npar);
  
  dpar = par / control.parscale;
  
  if (method_ == "Nelder-Mead") {
    nmmin(npar, dpar.memptr(), opar.memptr(), &val_, fminfn, &fail_,
          control.abstol, control.reltol, &func, control.alpha, control.beta,
          control.gamma, control.trace, &fncount_, control.maxit);
    
    par = opar % control.parscale;
    grcount_ = 0;
    
  } else if (method_ == "SANN") {
    int trace = control.trace;
    if (trace) trace = control.REPORT;
    
    if (control.tmax == NA_INTEGER || control.tmax < 1)
      Rcpp::stop("'tmax' is not a positive integer");
    
    internal::samin(npar, dpar.memptr(), &val_, fminfn, control.maxit,
                    control.tmax, control.temp, trace, &func);
    par = dpar % control.parscale;
    fncount_ = npar > 0 ? control.maxit : 1;
    grcount_ = 0;
    
  } else if (method_ == "BFGS") {
    arma::ivec mask = arma::ones<arma::ivec>(npar);
    vmmin(npar, dpar.memptr(), &val_, fminfn, fmingr, control.maxit,
          control.trace, mask.memptr(), control.abstol, control.reltol,
          control.REPORT, &func, &fncount_, &grcount_, &fail_);
    
    par = dpar % control.parscale;
  } else if (method_ == "CG") {
    cgmin(npar, dpar.memptr(), opar.memptr(), &val_, fminfn, fmingr, &fail_,
          control.abstol, control.reltol, &func, control.type, control.trace,
          &fncount_, &grcount_, control.maxit);
    
    par = opar % control.parscale;
  } else if (method_ == "L-BFGS-B") {
    arma::vec lower(npar);
    arma::vec upper(npar);
    arma::ivec nbd = arma::zeros<arma::ivec>(npar);
    char msg[60];
    
    for (std::size_t i = 0; i != npar; ++i) {
      lower(i) = lower_(i) / func.os.parscale_(i);
      upper(i) = upper_(i) / func.os.parscale_(i);
      if (!std::isfinite(lower(i))) {
        if (!std::isfinite(upper(i)))
          nbd(i) = 0;
        else
          nbd(i) = 3;
      } else {
        if (!std::isfinite(upper(i)))
          nbd(i) = 1;
        else
          nbd(i) = 2;
      }
    }
    
    func.os.usebounds_ = 1;
    func.os.lower_ = lower;
    func.os.upper_ = upper;
    
    lbfgsb(npar, control.lmm, dpar.memptr(), lower.memptr(), upper.memptr(),
           nbd.memptr(), &val_, fminfn, fmingr, &fail_, &func, control.factr,
           control.pgtol, &fncount_, &grcount_, control.maxit, msg,
           control.trace, control.REPORT);
    
    par = dpar % control.parscale;
    message_ = msg;
  } else
    Rcpp::stop("Roptim::minimize(): unknown 'method'");
  
  par_ = par;
  val_ *= func.os.fnscale_;
  
  if (hessian_flag_) func.ApproximateHessian(par_, hessian_);
}

}  // namespace roptim

#endif  // ROPTIM_H_


