//
//  loss.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OBJECTIVE_LOSS_HPP_
#define NSOPTIM_OBJECTIVE_LOSS_HPP_

namespace nsoptim {

//! Boilerplate base class for all loss functions.
//!
//! Loss functions must at least the following methods:
//! `data()` to give read access to the internal data, and
//! `operator(where)` to evaluate the loss at the given coefficients values.
//! `ZeroCoefficients()` to obtain the 0-coefficient value.
//!
//! Loss functions can optionally also implement the following methods:
//! `Difference(a, b)` to evaluate the difference of two coefficient values.
//!
//! Loss functions should be easy and quick to copy and move. The main purpose is not to provide functionality but
//! context.
template<class Data>
class LossFunction {
 public:
  using DataType = Data;

//! Access the data the loss operates on.
  //!
  //! @return the data the loss operates on.
  //! const Data& data() const;

  //! Evaluate the loss function.
  //!
  //! @param where where to evaluate the loss function.
  //! @return the loss evaluated at the given coefficients.
  //! double operator()(const Coefficients& where) const;
  //! Get the zero coefficients for this loss type.
  //!
  //! @return zero coefficients.
  // Coefficients ZeroCoefficients() const;

  //! Get the difference between two sets of coefficients.
  //!
  //! @param x a set of regression coefficients.
  //! @param y the other set of regression coefficients.
  //! @return the relative difference between `x` and `y`.
  // double Difference(const Coefficients& x, const Coefficients& y) const;
};
}  // namespace nsoptim

#endif  // NSOPTIM_OBJECTIVE_LOSS_HPP_
