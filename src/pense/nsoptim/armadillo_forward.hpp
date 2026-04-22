//
//  armadillo_forward.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_ARMADILLO_FORWARD_HPP_
#define NSOPTIM_ARMADILLO_FORWARD_HPP_

#define ARMA_DONT_USE_OPENMP 1

#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Weverything"
#endif

# include <RcppArmadillo/interface/RcppArmadilloForward.h>

# ifndef Rcpp__exceptions_impl__h
#  define Rcpp__exceptions_impl__h
namespace Rcpp {
inline void exception::record_stack_trace() {}
inline void exception::copy_stack_trace_to_r() const { rcpp_set_stack_trace(R_NilValue); }
}  // namespace Rcpp
# endif

# ifndef Rcpp_protection_meat_H
#  define Rcpp_protection_meat_H
namespace Rcpp {
template <typename T>
template <typename U>
Armor<T>::Armor(U x) : data() {
    init(wrap(x));
}

template <typename T>
template <typename U>
inline Armor<T>& Armor<T>::operator=(const U& x) {
    REPROTECT(data = wrap(x), index);
    return *this;
}
}  // namespace Rcpp
# endif

#ifdef __clang__
#  pragma clang diagnostic pop
#endif

#endif  // NSOPTIM_ARMADILLO_FORWARD_HPP_
