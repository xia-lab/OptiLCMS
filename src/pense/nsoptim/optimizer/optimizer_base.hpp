//
//  optimizer_base.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2019-01-02.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_OPTIMIZER_BASE_HPP_
#define NSOPTIM_OPTIMIZER_OPTIMIZER_BASE_HPP_

#include "optimum.hpp"
#include "../traits/traits.hpp"

namespace nsoptim {
//! Base class for all optimizer using loss function type `T`, penalty function type `U` and coefficient type `V`.
//! This class checks whether `T` is a valid loss function for coefficient type `V` as well as if `U` is a valid
//! penalty function for coefficient type `V`.
template<typename T, typename U, typename V>
class Optimizer {
 public:
  using LossFunction = T;     //< Loss function type
  using PenaltyFunction = U;  //< Penalty function type
  using Coefficients = V;     //< Coefficients type
  using Optimum = nsoptim::Optimum<LossFunction, PenaltyFunction, Coefficients>;

  static_assert(traits::is_loss_function<LossFunction>::value,
                "LossFunction does not implement the loss function interface");
  static_assert(traits::is_penalty_function<PenaltyFunction>::value,
                "PenaltyFunction does not implement the penalty function interface");
  static_assert(traits::loss_supports_evaluation<LossFunction, Coefficients>::value,
                "LossFunction does not support evaluation of the coefficients.");
  static_assert(traits::can_evaluate<PenaltyFunction, Coefficients>::value,
                "PenaltyFunction does not support evaluation of the coefficients.");
};
}  // nsoptim

#endif  // NSOPTIM_OPTIMIZER_OPTIMIZER_BASE_HPP_
