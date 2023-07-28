//
//  forward.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OBJECTIVE_FORWARD_HPP_
#define NSOPTIM_OBJECTIVE_FORWARD_HPP_

#include "../config.hpp"
#include "../armadillo_forward.hpp"

namespace nsoptim {

//! Full definition at nsoptim/objective/ls_regression_loss.hpp
class WeightedLsRegressionLoss;
class LsRegressionLoss;

//! Full definition at nsoptim/objective/en_penalty.hpp
class EnPenalty;
class LassoPenalty;
class RidgePenalty;

//! Full definition at nsoptim/objective/adaptive_en_penalty.hpp
class AdaptiveEnPenalty;
class AdaptiveLassoPenalty;

}  // namespace nsoptim

#endif  // NSOPTIM_OBJECTIVE_FORWARD_HPP_
