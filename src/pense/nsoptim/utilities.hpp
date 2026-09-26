//
//  utilities.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2019-01-02.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_UTILITIES_HPP_
#define NSOPTIM_UTILITIES_HPP_

#include <cstddef>
#include <ostream>

namespace nsoptim {
//! Globally unique ID for all objects used in nsoptim.
class ObjectId {
 public:
  ObjectId() noexcept : id_(ObjectId::NextId()) {}
  //! Copying is allowed and will preserve the ID!
  ObjectId(const ObjectId&) = default;
  ObjectId& operator=(const ObjectId&) = default;
  //! Moving is allowed and will preserve the ID!
  ObjectId(ObjectId&&) = default;
  ObjectId& operator=(ObjectId&&) = default;

  //! Compare two IDs.
  //!
  //! @param other the other ID.
  //! @return true if the IDs are not equal.
  bool operator!=(const ObjectId& other) const noexcept {
    return id_ != other.id_;
  }

  //! Compare two IDs.
  //!
  //! @param other the other ID.
  //! @return true if the IDs are equal.
  bool operator==(const ObjectId& other) const noexcept {
    return id_ == other.id_;
  }

  friend std::ostream& operator<< (std::ostream& stream, const ObjectId& id) {
    stream << "0x" << std::hex << id.id_ << std::dec;
    return stream;
  }

  static ObjectId null() noexcept {
    return ObjectId(kNullId);
  }
 private:
  static constexpr std::size_t kNullId = 0;
  std::size_t id_;

  ObjectId(const std::size_t id) noexcept : id_(id) {}

  static std::size_t NextId() noexcept {
    static std::size_t obj_counter = kNullId;

    std::size_t next_id = kNullId;
    #pragma omp atomic capture
    next_id = ++obj_counter;

    return next_id;
  }
};
} // namespace nsoptim

#endif  // NSOPTIM_UTILITIES_HPP_
