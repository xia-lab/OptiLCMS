//
//  regularization_path.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef REGULARIZATION_PATH_NEW_HPP_
#define REGULARIZATION_PATH_NEW_HPP_

#include <memory>
#include <tuple>
#include <type_traits>

#include "nsoptim.hpp"

#include "alias.hpp"
#include "m_loss.hpp"
#include "omp_utils.hpp"

namespace pense {
namespace regpath {

//! Test two coefficients for approximate equivalence.
//!
//! The function assumes that the coefficient vectors are of the same dimension.
//!
//! @param a coefficients to be compared
//! @param b coefficients to be compared
//! @param eps numerical tolerance for comparison
//! @return true if the two coefficient vectors are approximately equivalent, false otherwise.
template<class Coefficients>
bool CoefficientsEquivalent(const Coefficients& a, const Coefficients& b,
                            const double eps) noexcept {
  const double int_diff = a.intercept - b.intercept;
  if (int_diff * int_diff < a.beta.n_elem * eps) {
    // The intercept is similar. Check if the slope is also similar.
    const double beta_diff = arma::norm(a.beta - b.beta, 2);
    if (int_diff * int_diff + beta_diff * beta_diff < eps) {
      // The slope is also similar.
      return true;
    }
  }
  return false;
}

//! A list of starting points with associated optimizer.
template<class Ordering, typename... Ts>
class OrderedTuples {
 public:
  using Element = std::tuple<Ts...>;

  enum class InsertResult { kGood, kBad, kDuplicate };

  //! Create a new list of unique elements.
  explicit OrderedTuples() noexcept
    : max_size_(0), order_(), size_(0), elements_() {}

  //! Create a new list of unique items of unlimited size using the given ordering.
  //!
  //! @param order ordering of the elements.
  explicit OrderedTuples(const Ordering& order) noexcept
    : max_size_(0), order_(order), size_(0), elements_() {}

  //! Create a new list of unique items of limited size.
  //!
  //! @param max_size maximum number of coefficients retained.
  explicit OrderedTuples(const size_t max_size) noexcept
    : max_size_(max_size), order_(), size_(0), elements_() {}

  //! Create a new list of unique items of limited size.
  //!
  //! @param max_size maximum number of coefficients retained.
  //! @param order instance of the ordering class.
  OrderedTuples(const size_t max_size, const Ordering& order) noexcept
    : max_size_(max_size), order_(order), size_(0), elements_() {}

  //! Default copy constructor. Copy assignment is not possible.
  OrderedTuples(const OrderedTuples&) = default;
  OrderedTuples& operator=(const OrderedTuples&) = delete;

  //! Move constructor.
  OrderedTuples(OrderedTuples&& other) noexcept :
      max_size_(other.max_size_), order_(std::move(other.order_)),
      size_(other.size_), elements_(std::move(other.elements_)) {
    other.size_ = 0;
  }
  OrderedTuples& operator=(OrderedTuples&&) = delete;

  //! Get the number of elements in the list.
  size_t Size() const noexcept {
    return size_;
  }

  //! Clear all elements from the list.
  void Clear() noexcept {
    elements_.clear();
    size_ = 0;
  }

  //! Insert coefficients into the container.
  // template<typename... Args>
  InsertResult Emplace(Ts&&... args) {
    // Check if the optimum's objective function value is good enough.
    if (max_size_ > 0 && size_ >= max_size_ &&
          order_.after(elements_.front(), std::forward<Ts>(args)...)) {
      return InsertResult::kBad;
    }

    // Determine insert position.
    const auto elements_end = elements_.end();
    auto current_it = elements_.begin();
    auto insert_it = elements_.before_begin();

    while (current_it != elements_end) {
      const bool better = order_.before(*current_it, std::forward<Ts>(args)...);

      // Check if the coefficients are equal to the element at question.
      if (!better && !order_.after(*current_it, std::forward<Ts>(args)...) &&
          order_.equivalent(*current_it, std::forward<Ts>(args)...)) {
        return InsertResult::kDuplicate;
      } else if (!better) {
        break;
      }

      ++current_it;
      ++insert_it;
    }

    // No duplicate has been detected. Add after the insert position.
    elements_.emplace_after(insert_it, std::forward<Ts>(args)...);
    // Ensure that the size stays within the limits.
    if (++size_ > max_size_ && max_size_ > 0) {
      elements_.erase_after(elements_.before_begin());
      --size_;
    }
    return InsertResult::kGood;
  }

  const alias::FwdList<Element>& Elements() const noexcept {
    return elements_;
  }

  alias::FwdList<Element>& Elements() noexcept {
    return elements_;
  }

 private:
  const size_t max_size_;
  Ordering order_;
  size_t size_;
  alias::FwdList<Element> elements_;
};

template<class Coefficients>
class DuplicateCoefficients {
 public:
  explicit DuplicateCoefficients(const double eps) noexcept : eps_(eps) {}

  //! Does the existing element come before the new element?
  template<typename Element, typename... Args>
  bool before(const Element&, const Coefficients&, Args&&...) const noexcept {
    return false;
  }

  //! Does the existing element come before the new element?
  template<typename Element, typename... Args>
  bool after(const Element&, const Coefficients&, Args&&...) const noexcept {
    return false;
  }

  //! Is the existing element equivalent to new element?
  template<typename Element, typename... Args>
  bool equivalent(const Element& el, const Coefficients& coefs, Args&&...) const noexcept {
    return CoefficientsEquivalent(std::get<0>(el), coefs, eps_);
  }

 private:
  const double eps_;
};

template<class Optimizer>
class OptimaOrder {
  using Coefficients = typename Optimizer::Coefficients;
  using Optimum = typename Optimizer::Optimum;
 public:
  explicit OptimaOrder(const double eps) noexcept : eps_(eps) {}

  //! Does the existing element come before the new element?
  template<typename Element, typename... Args>
  bool before(const Element& el, const Coefficients& coefs, const double objf_value,
              Args&&...) const noexcept {
    return std::get<1>(el) > objf_value + eps_;
  }

  //! Does the existing element come before the new element?
  template<typename Element, typename... Args>
  bool before(const Element& el, const Optimum& opt, Args&&...) const noexcept {
    return std::get<0>(el).objf_value > opt.objf_value + eps_;
  }

  //! Does the existing element come before the new element?
  template<typename Element, typename... Args>
  bool after(const Element& el, const Coefficients& coefs, const double objf_value,
             Args&&...) const noexcept {
    return std::get<1>(el) < objf_value - eps_;
  }

  //! Does the existing element come before the new element?
  template<typename Element, typename... Args>
  bool after(const Element& el, const Optimum& opt, Args&&...) const noexcept {
    return std::get<0>(el).objf_value < opt.objf_value - eps_;
  }

  //! Is the existing element equivalent to new element?
  template<typename Element, typename... Args>
  bool equivalent(const Element& el, const Coefficients& coefs, Args&&...) const noexcept {
    return CoefficientsEquivalent(std::get<0>(el), coefs, eps_);
  }

  //! Is the existing element equivalent to new element?
  template<typename Element, typename... Args>
  bool equivalent(const Element& el, const Optimum& opt, Args&&...) const noexcept {
    return CoefficientsEquivalent(std::get<0>(el).coefs, opt.coefs, eps_);
  }

 private:
  const double eps_;
};

template<class Coefficients, typename... Ts>
using UniqueCoefficients = OrderedTuples<DuplicateCoefficients<Coefficients>, Coefficients, Ts...>;

template<class Optimizer, typename... Ts>
using UniqueStartPoints = OrderedTuples<OptimaOrder<Optimizer>, typename Optimizer::Coefficients,
                                        double, Optimizer, Ts...>;

template<class Optimizer, typename... Ts>
using UniqueOptima = OrderedTuples<OptimaOrder<Optimizer>, typename Optimizer::Optimum, Optimizer,
                                   Ts...>;

} // namespace regpath

template<class Optimizer>
class RegularizationPath {
  using LossFunction = typename Optimizer::LossFunction;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using Coefficients = typename Optimizer::Coefficients;
  using PenaltyList = alias::FwdList<PenaltyFunction>;
  using Optimum = typename Optimizer::Optimum;
  using IndividualStartingPoints = alias::FwdList<alias::FwdList<Coefficients>>;
  using UniqueCoefficients = regpath::UniqueCoefficients<Coefficients>;
  using UniqueCoefficientsOrder = regpath::DuplicateCoefficients<Coefficients>;
  using MetricsPtr = std::unique_ptr<nsoptim::Metrics>;
  using ExploredSolutions = regpath::UniqueStartPoints<Optimizer, MetricsPtr>;
  using ExploredSolutionsOrder = regpath::OptimaOrder<Optimizer>;
  using BestOptima = regpath::UniqueOptima<Optimizer>;
  using BestOptimaOrder = regpath::OptimaOrder<Optimizer>;

 public:
  struct Solutions {
    const PenaltyFunction& penalty;
    alias::Optima<Optimizer> optima;
  };

  //! Create a regularization paths using the given optimizer, loss function, and list of penalties.
  //!
  //! @param optimizer the optimizer to use. The loss function has to be set already.
  //! @param penalties a list of penalty functions.
  //! @param max_optima the maximum number of optima per penalty level.
  //! @param comparison_tol numeric tolerance for comparing two optima.
  //! @param num_threads number of threads to use.
  RegularizationPath(const Optimizer& optimizer,
                     const PenaltyList& penalties, const int max_optima,
                     const double comparison_tol, const int num_threads) :
      optimizer_template_(optimizer), penalties_(penalties),
      max_optima_(max_optima), comparison_tol_(comparison_tol), num_threads_(num_threads),
      shared_starts_(UniqueCoefficientsOrder(comparison_tol_)),
      best_starts_(max_optima, BestOptimaOrder(comparison_tol)),
      penalties_it_(penalties_.begin()) {
        auto penalties_it = penalties_.before_begin();
        const auto penalties_end = penalties_.end();

        while (++penalties_it != penalties_end) {
          individual_starts_.emplace_front(
            UniqueCoefficients(UniqueCoefficientsOrder(comparison_tol_)));
        }
        individual_starts_it_ = individual_starts_.before_begin();
      }

  //! Set the exploration options.
  //!
  //! @param explore_it the number of iterations for exploration. If <= 0, no exploration will
  //!   be done and all starting points will be iterated to full convergence.
  //! @param explore_tol the numeric tolerance for exploring solutions.
  //! @param explored_keep how many explored solutions to keep for full concentration.
  void ExplorationOptions(const int explore_it, const double explore_tol,
                          const int explored_keep) noexcept {
    explore_it_ = explore_it;
    explore_tol_ = explore_tol;
    explored_keep_ = explored_keep;
  }

  //! Enable/disable carrying forward solutions from the previous penalty.
  //!
  //! @param enabled whether to enable warm starts or not.
  void EnableWarmStarts(const bool enabled) noexcept {
    use_warm_start_ = enabled;
  }

  //! Add a starting point to be used only at the specified penalty.
  //!
  //! @param penalty penalty at which the starting point should be used.
  //! @param coefs starting point.
  void EmplaceIndividualStartingPoints(IndividualStartingPoints&& coefs_lists) {
    auto emplace_it = individual_starts_.begin();
    for (auto&& coefs_list : coefs_lists) {
      for (auto&& coefs : coefs_list) {
        emplace_it->Emplace(std::move(coefs));
      }
      emplace_it++;
    }
  }

  //! Add a starting point to be used for all penalties.
  //!
  //! @param coefs starting point.
  void EmplaceSharedStartingPoint(Coefficients&& coefs) {
    shared_starts_.Emplace(std::move(coefs));
  }

  Solutions Next() {
    ++individual_starts_it_;
    const auto& current_penalty = *penalties_it_++;
    optimizer_template_.penalty(current_penalty);

    auto explored_solutions = explore_it_ > 0 ? Explore() : SkipExploration();
    return Solutions { current_penalty, Concentrate(std::move(explored_solutions)) };
  }

  bool End() const noexcept {
    return penalties_it_ == penalties_.end();
  }

 private:
  Optimizer optimizer_template_;
  const PenaltyList& penalties_;
  const int max_optima_;
  const double comparison_tol_;
  int num_threads_;  //< OpenMP requires it to be an lvalue!
  bool use_warm_start_ = true;
  int explore_it_ = 0;
  double explore_tol_ = 0;
  int explored_keep_ = 1;

  alias::FwdList<UniqueCoefficients> individual_starts_;
  UniqueCoefficients shared_starts_;
  BestOptima best_starts_;

  typename alias::FwdList<UniqueCoefficients>::iterator individual_starts_it_;
  typename PenaltyList::const_iterator penalties_it_;

  ExploredSolutions Explore() {
    if (omp::Enabled(num_threads_)) {
      return Explore(std::true_type{});
    } else {
      return Explore(std::false_type{});
    }
  }

  ExploredSolutions Explore(std::true_type) {
    ExploredSolutions explored_solutions(explored_keep_, ExploredSolutionsOrder(comparison_tol_));
    const auto is_end = individual_starts_it_->Elements().end();
    const auto sh_end = shared_starts_.Elements().end();

    #pragma omp parallel \
                num_threads(num_threads_) \
                default(shared)
    {
      #pragma omp single nowait
      for (auto is_it = individual_starts_it_->Elements().begin(); is_it != is_end; ++is_it) {
        #pragma omp task \
                    default(none) \
                    firstprivate(is_it) \
                    shared(explore_tol_, explore_it_) \
                    shared(explored_solutions, optimizer_template_)
        {
          Optimizer optimizer(optimizer_template_);
          optimizer.convergence_tolerance(explore_tol_);
          auto optimum = optimizer.Optimize(std::get<0>(*is_it), explore_it_);

          #pragma omp critical(insert_explored)
          explored_solutions.Emplace(std::move(optimum.coefs), std::move(optimum.objf_value),
                                     std::move(optimizer), std::move(optimum.metrics));

        }
      }

      #pragma omp single nowait
      for (auto sh_it = shared_starts_.Elements().begin(); sh_it != sh_end; ++sh_it) {
        #pragma omp task \
                    firstprivate(sh_it) \
                    default(none) \
                    shared(explore_tol_, explore_it_) \
                    shared(explored_solutions, optimizer_template_)
        {
          Optimizer optimizer(optimizer_template_);
          optimizer.convergence_tolerance(explore_tol_);
          auto optimum = optimizer.Optimize(std::get<0>(*sh_it), explore_it_);

          #pragma omp critical(insert_explored)
          explored_solutions.Emplace(std::move(optimum.coefs), std::move(optimum.objf_value),
                                     std::move(optimizer), std::move(optimum.metrics));

        }
      }

      #pragma omp single nowait
      if (use_warm_start_ || explored_solutions.Size() == 0) {
      const auto bs_end = best_starts_.Elements().end();

        for (auto bs_it = best_starts_.Elements().begin(); bs_it != bs_end; ++bs_it) {
          #pragma omp task \
                      firstprivate(bs_it) \
                      default(none) \
                      shared(explore_tol_, explore_it_, explored_solutions, optimizer_template_)
          {
            auto&& optimizer = std::get<1>(*bs_it);
            optimizer.convergence_tolerance(explore_tol_);
            optimizer.penalty(optimizer_template_.penalty());
            auto optimum = optimizer.Optimize(explore_it_);

            #pragma omp critical(insert_explored)
            explored_solutions.Emplace(std::move(optimum.coefs), std::move(optimum.objf_value),
                                       std::move(optimizer), std::move(optimum.metrics));

          }
        }
      }
    }

    Rcpp::checkUserInterrupt();
    return explored_solutions;
  }

  ExploredSolutions Explore(std::false_type) {
    ExploredSolutions explored_solutions(explored_keep_, ExploredSolutionsOrder(comparison_tol_));

    for (auto& start : individual_starts_it_->Elements()) {
      Optimizer optimizer(optimizer_template_);
      optimizer.convergence_tolerance(explore_tol_);
      auto optimum = optimizer.Optimize(std::get<0>(start), explore_it_);
      explored_solutions.Emplace(std::move(optimum.coefs), std::move(optimum.objf_value),
                                 std::move(optimizer), std::move(optimum.metrics));

      Rcpp::checkUserInterrupt();
    }

    for (auto& start : shared_starts_.Elements()) {
      Optimizer optimizer(optimizer_template_);
      optimizer.convergence_tolerance(explore_tol_);
      auto optimum = optimizer.Optimize(std::get<0>(start), explore_it_);
      explored_solutions.Emplace(std::move(optimum.coefs), std::move(optimum.objf_value),
                                 std::move(optimizer), std::move(optimum.metrics));

      Rcpp::checkUserInterrupt();
    }

    if (use_warm_start_ || explored_solutions.Size() == 0) {
      for (auto& start : best_starts_.Elements()) {
        auto&& optimizer = std::get<1>(start);
        optimizer.convergence_tolerance(explore_tol_);
        optimizer.penalty(optimizer_template_.penalty());
        auto optimum = optimizer.Optimize(explore_it_);
        explored_solutions.Emplace(std::move(optimum.coefs), std::move(optimum.objf_value),
                                  std::move(optimizer), std::move(optimum.metrics));

        Rcpp::checkUserInterrupt();
      }
    }
    return explored_solutions;
  }

  //! Simply add all starts for the current penalty
  ExploredSolutions SkipExploration() {
    ExploredSolutions explored_solutions(0, ExploredSolutionsOrder(comparison_tol_));

    for (auto& start : individual_starts_it_->Elements()) {
      explored_solutions.Emplace(std::move(std::get<0>(start)), -1, Optimizer(optimizer_template_),
                                 MetricsPtr());
    }

    for (auto&& start : shared_starts_.Elements()) {
      // Shared starts must be copied, not moved!
      auto coefs = std::get<0>(start);
      explored_solutions.Emplace(std::move(coefs), -1, Optimizer(optimizer_template_),
                                 MetricsPtr());
    }

    if (use_warm_start_ || explored_solutions.Size() == 0) {
      for (auto& start : best_starts_.Elements()) {
        auto&& optimizer = std::get<1>(start);
        optimizer.penalty(optimizer_template_.penalty());
        explored_solutions.Emplace(std::move(std::get<0>(start).coefs), -1, std::move(optimizer),
                                  MetricsPtr());
      }
    }
    return explored_solutions;
  }

  alias::Optima<Optimizer> Concentrate(ExploredSolutions&& explored) {
    best_starts_.Clear();

    if (omp::Enabled(num_threads_)) {
      Concentrate(std::move(explored), std::true_type{});
    } else {
      Concentrate(std::move(explored), std::false_type{});
    }

    alias::Optima<Optimizer> optima;
    for (auto&& element : best_starts_.Elements()) {
      optima.emplace_front(std::get<0>(element));
    }
    return optima;
  }

  void Concentrate(ExploredSolutions&& explored, std::false_type) {
    const double conv_threshold = optimizer_template_.convergence_tolerance();

    for (auto&& start : explored.Elements()) {
      auto&& optimizer = std::get<2>(start);
      optimizer.convergence_tolerance(conv_threshold);
      auto optim = (std::get<1>(start) > 0) ?
        optimizer.Optimize() :
        optimizer.Optimize(std::get<0>(start));

      if (optim.metrics && std::get<3>(start)) {
        auto&& exploration_metrics = optim.metrics->CreateSubMetrics("exploration");
        exploration_metrics.AddSubMetrics(std::move(*std::get<3>(start)));
        std::get<3>(start).reset();
      }
      best_starts_.Emplace(std::move(optim), std::move(optimizer));

      Rcpp::checkUserInterrupt();
    }
  }

  void Concentrate(ExploredSolutions&& explored, std::true_type) {
    const double conv_threshold = optimizer_template_.convergence_tolerance();
    const auto ex_end = explored.Elements().end();

    #pragma omp parallel \
                num_threads(num_threads_) \
                default(shared)
    {
      #pragma omp single nowait
      for (auto ex_it = explored.Elements().begin(); ex_it != ex_end; ++ex_it) {
        #pragma omp task \
                    default(none) \
                    firstprivate(ex_it, conv_threshold) \
                    shared(best_starts_)
        {
          auto&& optimizer = std::get<2>(*ex_it);
          optimizer.convergence_tolerance(conv_threshold);
          auto optim = (std::get<1>(*ex_it) > 0) ?
            optimizer.Optimize() :
            optimizer.Optimize(std::get<0>(*ex_it));

          if (optim.metrics && std::get<3>(*ex_it)) {
            auto&& exploration_metrics = optim.metrics->CreateSubMetrics("exploration");
            exploration_metrics.AddSubMetrics(std::move(*std::get<3>(*ex_it)));
            std::get<3>(*ex_it).reset();
          }
          #pragma omp critical(insert_concentrated)
          best_starts_.Emplace(std::move(optim), std::move(optimizer));
        }
      }
    }
    Rcpp::checkUserInterrupt();
  }
};
} // namespace pense

#endif // REGULARIZATION_PATH_NEW_HPP_
