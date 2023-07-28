//
//  omp_utils.hpp
//  pense
//
//  Created by David Kepplinger on 2019-11-02.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef OMP_UTILS_HPP_
#define OMP_UTILS_HPP_

#include "autoconfig.hpp"

#ifdef PENSE_DISABLE_OPENMP
#  undef PENSE_ENABLE_OPENMP
#endif

#ifdef PENSE_ENABLE_OPENMP
#  undef PENSE_DISABLE_OPENMP
#endif

#define const_shared(list)

#ifdef PENSE_ENABLE_OPENMP

#include <omp.h>

namespace pense {
namespace omp {

#ifdef PENSE_OPENMP_ADD_CONST_SHARED
#   undef const_shared
#   define const_shared(list) shared(list)
#endif

//! Returns ``true`` if OpenMP is enabled.
inline bool Enabled(const int nr_threads) noexcept {
  return nr_threads > 1;
}

//! A conditional lock.
//! The lock is only active, if it is constructed as such.
class Lock {
 public:
  //! A lock which is only activated if the first argument is set to ``true``,
  inline explicit Lock(const bool enabled = true) noexcept : enabled_(enabled) {
    if (enabled_) {
      omp_init_lock(&lock_);
    }
  }

  //! A lock can not be copied, moved, or assigned to!
  Lock(const Lock&) = delete;
  Lock(Lock&&) = delete;
  Lock& operator=(const Lock&) = delete;
  Lock& operator=(Lock&&) = delete;

  virtual ~Lock() {
    if (enabled_) {
      omp_destroy_lock(&lock_);
    }
  }

  //! Acquire the lock.
  inline void Acquire() noexcept {
    if (enabled_) {
      omp_set_lock(&lock_);
    }
  }

  //! Release the lock.
  inline void Release() noexcept {
    if (enabled_) {
      omp_unset_lock(&lock_);
    }
  }

 private:
  const bool enabled_;
  omp_lock_t lock_;
};

//! An implicit guard which locks the given lock on construction and unlocks the lock when the guard goes out of scope.
class Guard {
 public:
  explicit Guard(Lock* lock) noexcept : lock_(lock) {
    lock_->Acquire();
  }

  virtual ~Guard() noexcept {
    lock_->Release();
  }

 private:
  Lock* lock_;
};

}  // namespace omp
}  // namespace pense

#else

namespace pense {
namespace omp {

//! Return ``true` if OpenMP is enabled.
inline constexpr bool Enabled(const int) noexcept {
  return false;
}

//! A lock object.
//! If OpenMP support is disabled, this is just a dummy which does not do anything.
class Lock {
 public:
  explicit Lock(const bool = true) noexcept {}
  void Acquire() const noexcept {}
  void Release() const noexcept {}
};

class Guard {
 public:
  explicit Guard(Lock* lock) noexcept {}
};

}  // namespace omp
}  // namespace pense

#endif

#endif  // OMP_UTILS_HPP_
