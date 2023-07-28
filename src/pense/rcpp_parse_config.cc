//
//  rcpp_parse_config.cc
//  pense
//
//  Created by David Kepplinger on 2019-05-01.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//
#include "nsoptim.hpp"

#include "rcpp_parse_config.hpp"
#include "rcpp_utils.hpp"
#include "cd_pense.hpp"

#include "nsoptim.hpp"

namespace {
constexpr int kAdmmMaxIt = 1000;
constexpr double kAdmmAcceleration = 1;

constexpr int kCDLsMaxIt = 1000;
constexpr int kCDLsResetIt = 8;

constexpr int kCDPenseMaxIt = 1000;
constexpr int kCDPenseResetIt = 8;
constexpr double kCDPenseLinesearchMult = 0.;
constexpr int kCDPenseLinesearchSteps = 10;

constexpr int kDalMaxIt = 100;
constexpr int kDalMaxInnerIt = 100;
constexpr double kDalEtaMult = 2;
constexpr double kDalEtaStartNumeratorCons = 0.01;
constexpr double kDalEtaStartNumeratorAggr = 1;
constexpr double kDalLambdaRelChangeAggr = 0.25;

constexpr int kMmMaxIt = 500;
constexpr nsoptim::MMConfiguration::TighteningType kMmTightening = nsoptim::MMConfiguration::TighteningType::kAdaptive;
constexpr int kMmTighteningSteps = 10;
}  // namespace

namespace Rcpp {
namespace traits {

nsoptim::AdmmLinearConfiguration Exporter<nsoptim::AdmmLinearConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  nsoptim::AdmmLinearConfiguration tmp = {
    pense::GetFallback(config_list, "max_it", kAdmmMaxIt),
    pense::GetFallback(config_list, "accelerate", kAdmmAcceleration)
  };
  return tmp;
}

nsoptim::DalEnConfiguration Exporter<nsoptim::DalEnConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  nsoptim::DalEnConfiguration tmp = {
      pense::GetFallback(config_list, "max_it", kDalMaxIt),
      pense::GetFallback(config_list, "max_inner_it", kDalMaxInnerIt),
      pense::GetFallback(config_list, "eta_start_numerator_conservative", kDalEtaStartNumeratorCons),
      pense::GetFallback(config_list, "eta_start_numerator_aggressive", kDalEtaStartNumeratorAggr),
      pense::GetFallback(config_list, "lambda_relchange_aggressive", kDalLambdaRelChangeAggr),
      pense::GetFallback(config_list, "eta_multiplier", kDalEtaMult)
  };
  return tmp;
}

pense::CDPenseConfiguration Exporter<pense::CDPenseConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  pense::CDPenseConfiguration tmp = {
      pense::GetFallback(config_list, "max_it", kCDPenseMaxIt),
      pense::GetFallback(config_list, "linesearch_mult", kCDPenseLinesearchMult),
      pense::GetFallback(config_list, "linesearch_steps", kCDPenseLinesearchSteps),
      pense::GetFallback(config_list, "reset_it", kCDPenseResetIt)
  };
  return tmp;
}

nsoptim::CDConfiguration Exporter<nsoptim::CDConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  nsoptim::CDConfiguration tmp = {
      pense::GetFallback(config_list, "max_it", kCDLsMaxIt),
      pense::GetFallback(config_list, "reset_it", kCDLsResetIt)
  };
  return tmp;
}

nsoptim::MMConfiguration Exporter<nsoptim::MMConfiguration>::get() const {
  const Rcpp::List config_list = as<const Rcpp::List>(r_obj_);
  nsoptim::MMConfiguration tmp = {
      pense::GetFallback(config_list, "max_it", kMmMaxIt),
      pense::GetFallback(config_list, "tightening", kMmTightening),
      pense::GetFallback(config_list, "tightening_steps", kMmTighteningSteps)
  };
  return tmp;
}

}  // namespace traits
}  // namespace Rcpp
