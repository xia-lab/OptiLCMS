//
//  rcpp_integration.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_RCPP_INTEGRATION_HPP_
#define NSOPTIM_RCPP_INTEGRATION_HPP_

#include "armadillo.hpp"

namespace Rcpp {
//! Wrap a sparse vector into a Matrix::sparseVector object
template <typename T>
SEXP wrap(const arma::SpCol<T>& svec) {
  const int  RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;

  // important: update internal state of SpMat object
  svec.sync();
  IntegerVector length = IntegerVector::create(svec.n_elem);

  // copy the data into R objects
  const Vector<RTYPE> values(svec.values, svec.values + svec.n_nonzero);
  IntegerVector rowind(svec.row_indices, svec.row_indices + svec.n_nonzero);

  // the sparseVector uses 1-based row indices.
  for (arma::uword i = 0; i < svec.n_nonzero; ++i) {
    rowind[i] += 1;
  }

  S4 r_sparse_vector("dsparseVector");
  r_sparse_vector.slot("length") = length;
  r_sparse_vector.slot("i") = rowind;
  r_sparse_vector.slot("x") = values;
  return r_sparse_vector;
}

//! Specialize Rcpp::wrap for nsoptim::RegressionCoefficients
template<typename T>
SEXP wrap(const nsoptim::RegressionCoefficients<T>& coefs) {
  return List::create(Named("intercept") = coefs.intercept,
                      Named("beta") = coefs.beta);
}

namespace traits {
//! Specialize Rcpp::as for nsoptim::RegressionCoefficients
template<typename T>
class Exporter<nsoptim::RegressionCoefficients<T>> {
 public:
  explicit Exporter(SEXP robj) {
    const List coef_list(robj);
    coefs_.intercept = Rcpp::as<double>(coef_list["intercept"]);
    coefs_.beta = Rcpp::as<T>(coef_list["beta"]);
  }

  nsoptim::RegressionCoefficients<T> get() const {
    return coefs_;
  }

 private:
  nsoptim::RegressionCoefficients<T> coefs_;
};

//! Specialize Rcpp::as for arma::SpCol<T>
template<typename T>
class Exporter<arma::SpCol<T>> {
 public:
  explicit Exporter(SEXP r_obj) {
    // Assume that the given R object is of type S4 (dsparseVector)
    S4 r_svec(r_obj);
    if (r_svec.is("dsparseVector")) {
      const auto nrows = as<arma::uword>(r_svec.slot("length"));
      const auto rowind = as<arma::uvec>(r_svec.slot("i"));
      SEXP val_slot = r_svec.slot("x");
      const arma::vec values(REAL(val_slot), Rf_length(val_slot), false, true);
      const arma::uvec colptr {0, rowind.n_elem};
      obj_ = arma::SpMat<T>(rowind - 1, colptr, values, nrows, 1).col(0);
    }
  }

  arma::SpMat<T> get() const {
    return obj_;
  }
 private:
  arma::SpMat<T> obj_;
};
}  // namespace traits
}  // namespace Rcpp

#endif  // NSOPTIM_RCPP_INTEGRATION_HPP_
