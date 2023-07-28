//
//  enpy_initest.cc
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include <algorithm>

#include "nsoptim.hpp"
#include "alias.hpp"
#include "constants.hpp"
#include "rcpp_utils.hpp"
#include "enpy_initest.hpp"

using arma::uvec;
using arma::mat;
using arma::uword;
using HashSet = std::unordered_set<uword>;
using subview_vec = arma::subview_col<double>;
using pense::alias::FwdList;

namespace {
constexpr int kDefaultCfgMaxIt = 1;  //!< Maximum number of iterations.
constexpr double kKeepPscProportion = 0.5;  //!< Proportion of observations to keep based on PSCs.
constexpr bool kUseResidualThreshold = false;  //!< Use a fixed threshold instead of a proportion
                                               //!< to screen observations based on their
                                               //!< residuals.
constexpr double kKeepResidualsProportion = 0.5;  //!< Proportion of observations to keep based
                                                  //!< on the residuals.
constexpr double kKeepResidualsThreshold = 2;  //!< Fixed threshold to keep observations based
                                               //!< on the residuals.
constexpr double kRetainBestFactor = 1.1;  //!< Retain not only the candidates from the last iteration,
                                           //!< but also those that are within this factor of the best candidate.
constexpr int kRetainMax = -1;  //!< Retain all candidates.
constexpr int kDefaultNumThreads = 1;  //!< Default number of threads.


inline uword HashUpdate(const uword hash, const uword value) noexcept;

//! Sort a range of values and compute the hash of the sorted range.
//! Because the range is first sorted, the hash is *independent* of the initial order.
//!
//! @param start iterator pointing to the start of the range.
//! @param end iterator pointing to the end of the range.
//! @return the hash of `vector`.
template<typename RandomAccessIterator>
uword SortAndHash(RandomAccessIterator start, RandomAccessIterator end) noexcept;

FwdList<uvec> GetSubsetList(const mat& pscs, const uvec& indices, const uword subset_size,
                            const bool use_indices) noexcept;

//! Comparator class to sort an _index_ vector based on ascending absolute values of the _values_ vector.
template<typename T>
class IndexCompAbsoluteAscending {
 public:
  //! Create a comparator using the given _values_ for comparision.
  explicit IndexCompAbsoluteAscending(const T& values) noexcept : values_(values) {}
  bool operator()(const uword a, const uword b) const {
    return std::abs(values_[a]) < std::abs(values_[b]);
  }
 private:
  const T& values_;
};

//! Comparator class to sort an _index_ vector based on ascending values of the _values_ vector.
template<typename T>
class IndexCompAscending {
 public:
  explicit IndexCompAscending(const T& values) noexcept : values_(values) {}
  bool operator()(const uword a, const uword b) const { return values_[a] < values_[b]; }
 private:
  const T& values_;
};

//! Comparator class to sort an _index_ vector based on descending values of the _values_ vector.
template<typename T>
class IndexCompDescending {
 public:
  explicit IndexCompDescending(const T& values) noexcept : values_(values) {}
  bool operator()(const uword a, const uword b) const { return values_[a] > values_[b]; }
 private:
  const T& values_;
};

}  // namespace

namespace pense {
namespace enpy_initest_internal {
PyConfiguration ParseConfiguration(const Rcpp::List& config) noexcept {
  return PyConfiguration{
    GetFallback(config, "max_it", kDefaultCfgMaxIt),
    GetFallback(config, "eps", kDefaultConvergenceTolerance),
    GetFallback(config, "keep_psc_proportion", kKeepPscProportion),
    GetFallback(config, "use_residual_threshold", kUseResidualThreshold),
    GetFallback(config, "keep_residuals_proportion", kKeepResidualsProportion),
    GetFallback(config, "keep_residuals_threshold", kKeepResidualsThreshold),
    GetFallback(config, "retain_best_factor", kRetainBestFactor),
    GetFallback(config, "retain_max", kRetainMax),
    GetFallback(config, "num_threads", kDefaultNumThreads)
  };
}

uword HashIndexVector(const uvec& vector) noexcept {
  uword hash = vector.n_elem;
  for (auto&& val : vector) {
    hash ^= HashUpdate(hash, val);
  }
  return hash;
}

uword HashSequence(const uword to) noexcept {
  uword hash = to + 1;
  for (uword val = 0; val <= to; ++val) {
    hash ^= HashUpdate(hash, val);
  }
  return hash;
}

arma::uvec GetResidualKeepIndices(const arma::vec& residuals, const double mscale_est,
                                  const PyConfiguration& config, arma::uvec* all_indices) {
  if (config.use_residual_threshold) {
    const double resid_threshold = config.keep_residuals_threshold * mscale_est;
    return arma::find(residuals <= resid_threshold);
  } else {
    const uword n_resid_keep = std::max<uword>(residuals.n_elem * config.keep_residuals_proportion,
                                                           kMinObs);
    std::partial_sort(all_indices->begin(), all_indices->begin() + n_resid_keep,
                      all_indices->end(), IndexCompAbsoluteAscending<arma::vec>(residuals));
    return arma::sort(all_indices->head_rows(n_resid_keep));
  }
}

alias::FwdList<uvec> GetSubsetList(const mat& pscs, const uword subset_size) noexcept {
  return ::GetSubsetList(pscs, uvec{}, subset_size, false);
}

alias::FwdList<uvec> GetSubsetList(const mat& pscs, const uvec& indices, const uword subset_size) noexcept {
  return ::GetSubsetList(pscs, indices, subset_size, true);
}
}  // namespace enpy_initest_internal
}  // namespace pense

namespace {
inline uword HashUpdate(const uword hash, const uword value) noexcept {
  return value + 0x9e3779b9 + (hash << 6) + (hash >> 2);
}

template<typename RandomAccessIterator>
uword SortAndHash(RandomAccessIterator start, RandomAccessIterator end) noexcept {
  uword hash = (end - start);
  std::sort(start, end);
  while (start != end) {
    hash ^= HashUpdate(hash, *start++);
  }
  return hash;
}

FwdList<uvec> GetSubsetList(const mat& pscs, const uvec& indices, const uword subset_size,
                            const bool use_indices) noexcept {
  FwdList<uvec> subsets;
  HashSet subset_candidate_hashes;
  uvec subset_indices = arma::regspace<uvec>(0, pscs.n_rows - 1);
  for (uword psc_col = 0; psc_col < pscs.n_cols; ++psc_col) {
    // Sort indices according to ascending absolute values of the PSC
    std::partial_sort(subset_indices.begin(), subset_indices.begin() + subset_size,
                      subset_indices.end(), IndexCompAbsoluteAscending<subview_vec>(pscs.col(psc_col)));

    const uword subset_abs_hash = SortAndHash(subset_indices.begin(), subset_indices.begin() + subset_size);
    if (subset_candidate_hashes.insert(subset_abs_hash).second) {
      // The subset is unique so add it to the list of subsets.
      if (use_indices) {
        subsets.emplace_front(indices.elem(subset_indices.head(subset_size)));
      } else {
        subsets.emplace_front(subset_indices.head(subset_size));
      }
    }

    // Sort indices according to ascending values of the PSC
    std::partial_sort(subset_indices.begin(), subset_indices.begin() + subset_size,
                      subset_indices.end(), IndexCompAscending<subview_vec>(pscs.col(psc_col)));

    const uword subset_asc_hash = SortAndHash(subset_indices.begin(), subset_indices.begin() + subset_size);
    if (subset_candidate_hashes.insert(subset_asc_hash).second) {
      // The subset is unique so add it to the list of subsets.
      if (use_indices) {
        subsets.emplace_front(indices.elem(subset_indices.head(subset_size)));
      } else {
        subsets.emplace_front(subset_indices.head(subset_size));
      }
    }

    // Sort indices according to descending values of the PSC
    std::partial_sort(subset_indices.begin(), subset_indices.begin() + subset_size,
                      subset_indices.end(), IndexCompDescending<subview_vec>(pscs.col(psc_col)));

    const uword subset_desc_hash = SortAndHash(subset_indices.begin(), subset_indices.begin() + subset_size);
    if (subset_candidate_hashes.insert(subset_desc_hash).second) {
      // The subset is unique so add it to the list of subsets.
      if (use_indices) {
        subsets.emplace_front(indices.elem(subset_indices.head(subset_size)));
      } else {
        subsets.emplace_front(subset_indices.head(subset_size));
      }
    }
  }

  return subsets;
}
}  // namespace
