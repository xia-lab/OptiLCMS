//
//  penalty.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OBJECTIVE_PENALTY_HPP_
#define NSOPTIM_OBJECTIVE_PENALTY_HPP_

namespace nsoptim {

//! Boilerplate base class for all penalty functions.
//!
//! Penalty functions must at least implement the following method:
//! `operator(where)` to evaluate the penalty at the given coefficients values.
//!
//! Penalty functions can optionally also implement the following methods:
//! `Difference(a, b)` to evaluate the difference of two coefficient values.
class PenaltyFunction {
 public:
  //! Evaluate the penalty function.
  //!
  //! @param where where to evaluate the penalty function.
  //! @return the penalty evaluated at the given coefficients.
  //! double operator()(const Coefficients& where) const;

  //! Get the difference between two sets of coefficients.
  //!
  //! @param x a set of regression coefficients.
  //! @param y the other set of regression coefficients.
  //! @return the relative difference between `x` and `y`.
  // double Difference(const Coefficients& x, const Coefficients& y) const;
};
}  // namespace nsoptim

#endif  // NSOPTIM_OBJECTIVE_PENALTY_HPP_
