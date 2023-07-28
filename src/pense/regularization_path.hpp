//
//  regularization_path.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef REGULARIZATION_PATH_HPP_
#define REGULARIZATION_PATH_HPP_

#include <memory>
#include <tuple>
#include <type_traits>

#include "nsoptim.hpp"

#include "alias.hpp"
#include "m_loss.hpp"
#include "omp_utils.hpp"

namespace pense {
namespace regularization_path {
template<class Optimizer, typename... Ts>
class UniqueOptima {
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using Optimum = typename Optimizer::Optimum;
  using Element = std::tuple<Optimum, Ts...>;

 public:
  enum class InsertResult { kGood, kBad, kDuplicate };

  UniqueOptima(const size_t max_size, const double eps) noexcept
      : max_size_(max_size), eps_(eps), size_(0), optima_() {}

  //! Default copy constructor. Copy assignment is not possible.
  UniqueOptima(const UniqueOptima&) = default;
  UniqueOptima& operator=(const UniqueOptima&) = delete;

  //! Move constructor.
  UniqueOptima(UniqueOptima&& other) noexcept : max_size_(other.max_size_), eps_(other.eps_), size_(other.size_),
                                                optima_(std::move(other.optima_)) {
    other.size_ = 0;
  }
  UniqueOptima& operator=(UniqueOptima&&) = delete;

  //! Check if an optimum is good enough to be potentially inserted into the container.
  bool IsGoodEnough(const Optimum& optimum) {
    // Check if the optimum's objective function value is good enough.
    if (size_ >= max_size_ && optimum.objf_value > std::get<0>(optima_.front()).objf_value) {
      return false;
    }
    return true;
  }

  //! Insert an optimum and associated information into the container, if the optimum is good enough.
  template<typename T, typename... Args>
  InsertResult Insert(T&& optimum, Args&&... args) {
    // Ensure that `T` is of type `Optimum`.
    static_assert(std::is_same<typename std::decay<T>::type, Optimum>::value, "Optimum is of wrong type.");

    // Check if the optimum's objective function value is good enough.
    if (size_ == max_size_ && optimum.objf_value > std::get<0>(optima_.front()).objf_value) {
      return InsertResult::kBad;
    }

    // Determine insert position.
    const auto optima_end = optima_.end();
    const auto optima_before_begin = optima_.before_begin();
    auto after_it = optima_.begin();
    auto insert_it = optima_before_begin;
    while (after_it != optima_end) {
      // Check if the two optima are equal.
      if (Equal(optimum, std::get<0>(*after_it))) {
        // Don't add.
        return InsertResult::kDuplicate;
      }

      // Check that the objective value of the given optimum is larger than the next optimum but smaller than
      // the current optimum.
      if ((optimum.objf_value > std::get<0>(*after_it).objf_value) &&
          (insert_it == optima_before_begin || std::get<0>(*insert_it).objf_value > optimum.objf_value)) {
        // Insert here.
        optima_.emplace_after(insert_it, std::forward<T>(optimum), std::forward<Args>(args)...);

        // Ensure that the size keeps within the limits.
        if (++size_ > max_size_) {
          optima_.erase_after(optima_.before_begin());
          --size_;
        }
        return InsertResult::kGood;
      }
      ++after_it;
      ++insert_it;
    }

    // All other optima have larger objective function value than the given optimum. Add at the end.
    optima_.emplace_after(insert_it, std::forward<T>(optimum), std::forward<Args>(args)...);
    // Ensure that the size keeps within the limits.
    if (++size_ > max_size_) {
      optima_.erase_after(optima_.before_begin());
      --size_;
    }
    return InsertResult::kGood;
  }

  const alias::FwdList<Element>& Elements() const noexcept {
    return optima_;
  }

  alias::FwdList<Element>& Elements() noexcept {
    return optima_;
  }

  template<size_t I>
  alias::FwdList<typename std::tuple_element<I, Element>::type > Elements() const noexcept {
    alias::FwdList<typename std::tuple_element<I, Element>::type > out;

    for (auto&& el : optima_) {
      out.push_front(std::get<I>(el));
    }

    return out;
  }

 private:
  const size_t max_size_;
  const double eps_;
  size_t size_;
  alias::FwdList<Element> optima_;

  bool Equal(const Optimum& a, const Optimum& b) const noexcept {
    // First check if the value of the objective function is similar.
    if (std::abs(a.objf_value - b.objf_value) < eps_) {
      // Value of the objective function is similar. Check if the intercept is similar.
      const double int_diff = a.coefs.intercept - b.coefs.intercept;
      if (int_diff * int_diff < eps_) {
        // The intercept is similar. Check if the slope is also similar.
        const double beta_diff = arma::norm(a.coefs.beta - b.coefs.beta, 2);
        if (beta_diff * beta_diff < eps_) {
          // The slope is also similar. Return true.
          return true;
        }
      }
    }
    return false;
  }
};

template<typename Optimizer>
using StartCoefficientsList = alias::FwdList<alias::FwdList<typename Optimizer::Coefficients>>;

//! A fixed-length list keeping track of optimizers and the value of the most recent optimum.
template<typename T>
class OptimizerList {
 private:
  using Optimum = typename T::Optimum;
  using Coefficients = typename T::Coefficients;
  using UOptima = UniqueOptima<T, T>;

 public:
  //! Initialize an optimizer-list with size `max_size` and numerical tolerance `eps`.
  //!
  //! @param max_size the maximum number of optimizers in the list.
  //! @param compare_tol_ the numerical tolerance for comparing two optima.
  OptimizerList(const int max_size, const double explore_tol, const int explore_it, const double compare_tol,
                const int num_threads) noexcept
      : max_size_(max_size), explore_tol_(explore_tol), explore_it_(explore_it), compare_tol_(compare_tol),
        num_threads_(num_threads), items_(max_size_, compare_tol) {}

  //! Add optimizer to the list.
  //!
  //! The list only retains the `max_size` best optima!
  //!
  //! @param optimizer the optimizer to compute the optimum starting from the different start points in `starts`.
  //! @param starts a list of start coefficients.
  //! @return a list of optima that have been attained from the given start points.
  alias::Optima<T> AddNew(T* optimizer, const alias::FwdList<Coefficients>& starts) {
    if (omp::Enabled(num_threads_)) {
      return AddNew(optimizer, starts, std::true_type{});
    } else {
      return AddNew(optimizer, starts, std::false_type{});
    }
  }

  //! Compute the unique optimum at the updated penalty.
  //!
  //! @param penalty the updated penalty function.
  //! @return a list of unique optima attained by the optimizers in the list.
  alias::Optima<T> UpdateAll(const typename T::PenaltyFunction& penalty) {
    if (omp::Enabled(num_threads_)) {
      UpdateAll(penalty, std::true_type{});
    } else {
      UpdateAll(penalty, std::false_type{});
    }
    return items_.template Elements<0>();
  }

 private:
  const int max_size_;
  const double explore_tol_;
  const int explore_it_;
  const double compare_tol_;
  int num_threads_;
  UOptima items_;

  //! UpdateAll if OpenMP support is enabled and needed.
  void UpdateAll(const typename T::PenaltyFunction& penalty, std::true_type) {
    UOptima old_optima = std::move(items_);
    #pragma omp parallel num_threads(num_threads_) shared(old_optima, items_, penalty) default(none)
    {
      #pragma omp single nowait
      for (auto item_it = old_optima.Elements().begin(), item_end = old_optima.Elements().end();
           item_it != item_end; ++item_it) {
        #pragma omp task firstprivate(item_it) default(shared)
        {
          auto&& optimizer = std::get<1>(std::move(*item_it));
          optimizer.penalty(penalty);
          auto&& optimum = optimizer.Optimize();

          #pragma omp critical(regpath_insert_optimum)
          items_.Insert(std::move(optimum), std::move(optimizer));
        }
      }
    }
  }

  //! UpdateAll if OpenMP support is disabled or not needed.
  void UpdateAll(const typename T::PenaltyFunction& penalty, std::false_type) {
    UOptima old_optima = std::move(items_);
    for (auto&& item : old_optima.Elements()) {
      auto&& optimizer = std::get<1>(item);
      optimizer.penalty(penalty);
      items_.Insert(optimizer.Optimize(), std::move(optimizer));
    }
  }

  //! AddNew if OpenMP support is enabled and needed.
  alias::Optima<T> AddNew(T* optimizer, const alias::FwdList<Coefficients>& starts, std::true_type) {
    // First briefly explore all starting points
    UOptima cold_items(max_size_, compare_tol_);
    alias::Optima<T> cold_optima;
    // Change to "approximating".
    const double original_tol = optimizer->convergence_tolerance();
    optimizer->convergence_tolerance(explore_tol_);

    #pragma omp parallel num_threads(num_threads_) default(none) shared(optimizer, starts, cold_items, cold_optima) \
      firstprivate(original_tol)
    {
      #pragma omp single
      for (auto start_it = starts.begin(), start_end = starts.end(); start_it != start_end; ++start_it) {
        #pragma omp task default(none) firstprivate(start_it) shared(optimizer, cold_items)
        {
          auto tmp_optimizer = *optimizer;
          auto cold_optimum = tmp_optimizer.Optimize(*start_it, explore_it_);
          if (cold_optimum.status != nsoptim::OptimumStatus::kError) {
            #pragma omp critical(regpath_insert_optimum)
            cold_items.Insert(std::move(cold_optimum), std::move(tmp_optimizer));
          }
        }
      }

      #pragma omp single nowait
      {
        // Fully iterate the best cold candidates and retain those that are better than the updated optima.
        for (auto cand_it = cold_items.Elements().begin(), cand_end = cold_items.Elements().end();
            cand_it != cand_end; ++cand_it) {
          #pragma omp task default(none) shared(cold_optima) firstprivate(cand_it, original_tol)
          {
            auto&& optimizer = std::get<1>(*cand_it);
            optimizer.convergence_tolerance(original_tol);
            auto&& optimum = optimizer.Optimize();

            #pragma omp critical(regpath_insert_optimum)
            {
              cold_optima.emplace_front(std::move(optimum));
              items_.Insert(cold_optima.front(), std::move(optimizer));
            }
          }
        }
      }
    }

    // Revert to original tolerance level.
    optimizer->convergence_tolerance(original_tol);

    return cold_optima;
  }

  //! AddNew if OpenMP support is disabled or not needed.
  alias::Optima<T> AddNew(T* optimizer, const alias::FwdList<Coefficients>& starts, std::false_type) {
    // First briefly explore all starting points
    UOptima cold_items(max_size_, compare_tol_);
    // Change to "approximating".
    const double original_tol = optimizer->convergence_tolerance();
    optimizer->convergence_tolerance(explore_tol_);

    for (auto&& start : starts) {
      auto tmp_optimizer = *optimizer;
      auto cold_optimum = tmp_optimizer.Optimize(start, explore_it_);
      if (cold_optimum.status != nsoptim::OptimumStatus::kError) {
        cold_items.Insert(std::move(cold_optimum), std::move(tmp_optimizer));
      }
    }

    // Revert to original tolerance level.
    optimizer->convergence_tolerance(original_tol);

    // Fully iterate the best cold candidates and retain those that are better than the updated optima.
    alias::Optima<T> cold_optima;
    for (auto&& candidate : cold_items.Elements()) {
      auto&& optimizer = std::get<1>(candidate);
      optimizer.convergence_tolerance(original_tol);
      cold_optima.emplace_front(optimizer.Optimize());
      const Optimum& cold_optimum = cold_optima.front();
      items_.Insert(cold_optimum, std::move(optimizer));
    }
    return cold_optima;
  }
};
}  // namespace regularization_path

//! 0-based Regularization Path
//!
//! Compute the regularization path for the same loss function but different penalties, starting the first optimization
//! at the 0-vector.
//! Subsequent solutions are based on starting the optimization at the previous solution.
template<typename Optimizer>
class RegPath0 {
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using PenaltyList = alias::FwdList<PenaltyFunction>;
  using Optimum = typename Optimizer::Optimum;

 public:
  //! Create the regularization path using the given optimizer, loss function, and list of penalties.
  //! @param optimizer the optimizer to use.
  //! @param loss the loss function to optimize.
  //! @param penalties a list of penalty functions.
  RegPath0(const Optimizer& optimizer, const LossFunction& loss, const PenaltyList& penalties) noexcept
      : penalties_(penalties), optim_(optimizer), penalty_it_(penalties.cbegin()) {
    optim_.loss(loss);
  }

  Optimum Next() {
    // At the first penalty, compute the solution from 0.
    if (penalty_it_ == penalties_.cbegin()) {
      const auto zero_coef = optim_.loss().template ZeroCoefficients<Coefficients>();
      optim_.penalty(*penalty_it_++);
      return optim_.Optimize(zero_coef);
    }
    // Otherwise, use the previous solution
    optim_.penalty(*penalty_it_++);
    return optim_.Optimize();
  }

  bool End() const noexcept {
    return penalty_it_ == penalties_.cend();
  }

 private:
  const PenaltyList& penalties_;
  Optimizer optim_;
  typename PenaltyList::const_iterator penalty_it_;
};

//! Regularization Path starting at the same starting point for every penalty.
//!
//! Compute the regularization path for the same loss function but different penalties.
//! At each penalty, the optimization is started from the same starting point!
template<typename Optimizer>
class RegPathIdentical {
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using PenaltyList = alias::FwdList<PenaltyFunction>;
  using Optimum = typename Optimizer::Optimum;

 public:
  //! Create the regularization path using the given optimizer, loss function, and list of penalties.
  //! @param optimizer the optimizer to use.
  //! @param loss the loss function to optimize.
  //! @param penalties a list of penalty functions.
  //! @param start the coefficients to start the optimization for every penalty.
  RegPathIdentical(const Optimizer& optimizer, const LossFunction& loss, const PenaltyList& penalties,
                   const Coefficients& start) noexcept
      : penalties_(penalties), start_(start), optim_(optimizer), penalty_it_(penalties.cbegin()), explored_(false) {
    optim_.loss(loss);
  }

  //! Compute the "exact" optimum at the current penalty level, then move on to the next penalty level.
  //! @return an "exact" optimum.
  Optimum Next() {
    if (explored_) {
      explored_ = false;
      // Penalty was incremented before!
      return optim_.Optimize();
    }
    // Increment only now!
    optim_.penalty(*penalty_it_++);
    return optim_.Optimize(start_);
  }

  //! Compute an approximate optimum using *larger than usual* tolerance at the current penalty level.
  //!
  //! @param eps the relaxed convergence tolerance.
  //! @return an approximate optimum.
  Optimum Explore(const double eps, const int maxit) {
    const double original_eps = optim_.convergence_tolerance();
    optim_.convergence_tolerance(eps);
    optim_.penalty(*penalty_it_);
    const Optimum tmp = optim_.Optimize(start_, maxit);
    optim_.convergence_tolerance(original_eps);

    // Increment
    ++penalty_it_;
    explored_ = true;

    return tmp;
  }

  //! Check if there are more penalties to compute the optimum for.
  bool End() const noexcept {
    return penalty_it_ == penalties_.cend();
  }

 private:
  const PenaltyList& penalties_;
  Coefficients start_;
  Optimizer optim_;
  typename PenaltyList::const_iterator penalty_it_;
  bool explored_;
};

//! Parallel Regularization Paths, starting at the best optima from the previous penalty.
//!
//! Compute the regularization path for the same loss function but different penalties.
//! At each penalty, several optimizations are performed:
//!  * starting from all solutions at the previous penalty (if available)
//!  * starting at all given starting points for the given penalty
//! From these solutions, only the best are retained and used in the subsequent optimizations.
template<typename Optimizer>
class RegPathCarryForward {
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using PenaltyList = alias::FwdList<PenaltyFunction>;
  using Optimum = typename Optimizer::Optimum;
  using StartCoefficients = alias::FwdList<alias::FwdList<Coefficients>>;

 public:
  class OptimumIterator;
  using iterator = OptimumIterator;

  //! Create the regularization path using the given optimizer, loss function, and list of penalties.
  //! @param loss the loss for which to compute the regularization path.
  //! @param penalties a list of penalty functions.
  //! @param optimizer the optimizer to use.
  //! @param starts a list with of starting points, i.e., for every item in `penalties`, `starts` contains a list of
  //!               coefficients. If the list is empty, the optimization starts at the best solutions at the previous
  //!               penalty.
  //! @param explore_tol convergence tolerance to approximate, or "explore", possible optima.
  //! @param explore_it maximum number of iterations to perform in the exploration state.
  //! @param nr_retain number of optima to keep track of.
  //! @param comparison_tol numerical tolerance to distinguish optima.
  RegPathCarryForward(const Optimizer& optimizer, const LossFunction& loss, const PenaltyList& penalties,
                      const StartCoefficients& starts, const double explore_tol, const int explore_it,
                      const int nr_retain, const double comparison_tol, const int num_threads) noexcept
    : penalties_(penalties), starts_(starts), optim_(optimizer), penalty_it_(penalties.cbegin()),
      starts_it_(starts_.cbegin()), skip_(starts_it_->empty()),
      cold_list_(nr_retain, explore_tol, explore_it, comparison_tol, num_threads) {
    optim_.loss(loss);
  }

  alias::Optima<Optimizer> Next() {
    if (skip_ && starts_it_->empty()) {
      ++penalty_it_;
      ++starts_it_;
      return alias::Optima<Optimizer>();
    }
    skip_ = false;

    // First update all previous solutions.
    auto optima = cold_list_.UpdateAll(*penalty_it_);

    // Now add new cold candidates.
    if (!starts_it_->empty()) {
      // Fully iterate the "best" optima.
      optim_.penalty(*penalty_it_);
      auto cold_optima = cold_list_.AddNew(&optim_, *starts_it_);
      optima.splice_after(optima.before_begin(), cold_optima);
    }
    ++penalty_it_;
    ++starts_it_;

    return optima;
  }

  bool End() const noexcept {
    return penalty_it_ == penalties_.cend();
  }

 private:
  const PenaltyList& penalties_;
  StartCoefficients starts_;
  Optimizer optim_;
  typename PenaltyList::const_iterator penalty_it_;
  typename StartCoefficients::const_iterator starts_it_;
  bool skip_;
  regularization_path::OptimizerList<Optimizer> cold_list_;
};

template<typename Optimizer>
class RegPathCombined {
 private:
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using PenaltyList = alias::FwdList<PenaltyFunction>;
  using Optimum = typename Optimizer::Optimum;
  using StartCoefficients = alias::FwdList<alias::FwdList<Coefficients>>;
  using UniqueOptima = regularization_path::UniqueOptima<Optimizer>;

  template<typename... Ts>
  using GenericUniqueOptima = regularization_path::UniqueOptima<Optimizer, Ts...>;

 public:
  //! Create a combination of regularization paths using the given optimizer, loss function, and list of penalties.
  //! @param optimizer the optimizer to use.
  //! @param loss the loss function to optimize.
  //! @param penalties a list of penalty functions.
  //! @param max_optima the maximum number of optima per penalty level.
  //! @param explore_tol the numeric tolerance for exploring solutions.
  //! @param comparison_tol numeric tolerance for comparing two optima.
  RegPathCombined(const Optimizer& optimizer, const LossFunction& loss, const PenaltyList& penalties,
                  const int max_optima, const int nr_explore, const double explore_tol, const int explore_it,
                  const double comparison_tol, const int num_threads) noexcept
      : optimizer_(optimizer), loss_(loss), penalties_(penalties), comparison_tol_(comparison_tol),
        max_optima_(max_optima), explore_tol_(explore_tol), explore_it_(explore_it), nr_explore_(nr_explore),
        num_threads_(num_threads) {}

  //! Add a 0-based regularization path for exploration.
  void Add() {
    if (!rp_0_) {
      rp_0_.reset(new RegPath0<Optimizer>(optimizer_, loss_, penalties_));
    }
  }

  //! Add a regularization path that always uses the given start coefficients, for all penalties.
  //! @param start the coefficients to start the optimization for every penalty.
  void Add(const Coefficients& start) {
    rp_id_.emplace_front(optimizer_, loss_, penalties_, start);
  }

  //! Add a regularization path that carries forward only the best optima.
  //! @param starts a list with of starting points, i.e., for every item in `penalties`, `starts` contains a list of
  //!               coefficients. If the list is empty, the optimization starts at the best solutions at the previous
  //!               penalty.
  //! @param nr_retain number of optima to keep track of.
  //! @param eps numerical tolerance to distinguish optima.
  void Add(const StartCoefficients& starts) {
    rp_cf_.emplace_front(optimizer_, loss_, penalties_, starts, explore_tol_, explore_it_, nr_explore_, comparison_tol_,
                         num_threads_);
  }

  //! Get the unique optima from all regularization paths at the next penalty level.
  alias::Optima<Optimizer> Next() {
    UniqueOptima next_optima(max_optima_, comparison_tol_);

    if (rp_0_) {
      next_optima.Insert(rp_0_->Next());
      Rcpp::checkUserInterrupt();
    }

    if (omp::Enabled(num_threads_)) {
      NextIdentical(&next_optima, std::true_type{});
    } else {
      NextIdentical(&next_optima, std::false_type{});
    }

    Rcpp::checkUserInterrupt();

    for (auto&& reg_path : rp_cf_) {
      for (auto&& optimum : reg_path.Next()) {
        next_optima.Insert(std::move(optimum));
        Rcpp::checkUserInterrupt();
      }
    }

    return next_optima.template Elements<0>();
  }

  //! Check if the regularization paths are at their end.
  bool End() const noexcept {
    if (rp_0_) {
      return rp_0_->End();
    }

    for (auto&& reg_path : rp_id_) {
      return reg_path.End();
    }

    for (auto&& reg_path : rp_cf_) {
      return reg_path.End();
    }

    return true;
  }

 private:
  const Optimizer& optimizer_;
  const LossFunction& loss_;
  const PenaltyList& penalties_;
  const double comparison_tol_;
  const int max_optima_;
  const double explore_tol_;
  const int explore_it_;
  const int nr_explore_;
  int num_threads_;  //!< Can not be constant, as OpenMP requires it to be an lvalue!
  std::unique_ptr<RegPath0<Optimizer>> rp_0_;
  alias::FwdList<RegPathIdentical<Optimizer>> rp_id_;
  alias::FwdList<RegPathCarryForward<Optimizer>> rp_cf_;

  //! Compute the next "identical" solutions if OpenMP support is enabled and needed.
  void NextIdentical(UniqueOptima* next_optima, std::true_type) {
    GenericUniqueOptima< RegPathIdentical<Optimizer>* > ident_explore_optima(nr_explore_, comparison_tol_);
    #pragma omp parallel num_threads(num_threads_) shared(rp_id_, next_optima, ident_explore_optima) default(none)
    {
      // First, explore all "identical" reg. paths.
      #pragma omp single
      for (auto reg_path_it = rp_id_.begin(), reg_path_end = rp_id_.end(); reg_path_it != reg_path_end; ++reg_path_it) {
        #pragma omp task firstprivate(reg_path_it) shared(ident_explore_optima) \
                         const_shared(explore_tol_), const_shared(explore_it_) default(none)
        {
          // This takes a while...
          auto&& optimum = reg_path_it->Explore(explore_tol_, explore_it_);
          if (ident_explore_optima.IsGoodEnough(optimum)) {
            // Try to add to the list of unique optima
            #pragma omp critical(insert_explored_optimum)
            ident_explore_optima.Insert(std::move(optimum), &(*reg_path_it));
          }
        }
      }

      // Then fully iterate the promising solutions.
      #pragma omp single nowait
      for (auto optimum_tuple_it = ident_explore_optima.Elements().begin(),
                optimum_tuple_end = ident_explore_optima.Elements().end(); optimum_tuple_it != optimum_tuple_end;
                ++optimum_tuple_it) {
        #pragma omp task firstprivate(optimum_tuple_it) shared(next_optima) default(none)
        {
          // This takes a while...
          auto&& optimum = std::get<1>(*optimum_tuple_it)->Next();
          // Add to the list of unique optima
          #pragma omp critical(insert_next_optima)
          next_optima->Insert(std::move(optimum));
        }
      }
    }
  }

  //! Compute the next "identical" solutions if OpenMP support is disabled or not needed.
  void NextIdentical(UniqueOptima* next_optima, std::false_type) {
    GenericUniqueOptima< RegPathIdentical<Optimizer>* > ident_explore_optima(nr_explore_, comparison_tol_);
    for (auto&& reg_path : rp_id_) {
      ident_explore_optima.Insert(reg_path.Explore(explore_tol_, explore_it_), &reg_path);
    }

    for (auto&& optimum_tuple : ident_explore_optima.Elements()) {
      next_optima->Insert(std::get<1>(optimum_tuple)->Next());
    }
  }
};
}  // namespace pense

#endif  // REGULARIZATION_PATH_HPP_
