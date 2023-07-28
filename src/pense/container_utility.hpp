//
//  container_utility.hpp
//  pense
//
//  Created by David Kepplinger on 2019-11-04.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef CONTAINER_UTILITY_HPP_
#define CONTAINER_UTILITY_HPP_

#include <functional>
#include <forward_list>

namespace pense {
namespace utility {
//! A std::forward_list with items ordered according to the value of another element.
template<typename T1, typename T2, typename Comparator = std::less<T1>>
class OrderedList {
  using ListType = std::forward_list<T2>;

 public:
  //! Create an empty ordered list.
  OrderedList() noexcept {}

  //! Create an empty ordered list using *comp* for comparisons.
  explicit OrderedList(const Comparator& comp) noexcept : comp_(comp) {}

  //! Insert an item at the position given by *order_item*.
  //!
  //! @return Iterator pointing to the inserted element.
  typename ListType::iterator insert(const T1& order_item, const T2& item) {
    auto order_it = order_items_.begin();
    auto insert_order_it = order_items_.before_begin();
    const auto order_end = order_items_.end();
    auto insert_item_it = items_.before_begin();

    while (order_it != order_end && comp_(*order_it, order_item)) {
      ++insert_item_it;
      ++insert_order_it;
      ++order_it;
    }
    order_items_.insert_after(insert_order_it, order_item);
    return items_.insert_after(insert_item_it, item);
  }

  //! Emplace an item at the position given by *order_item*.
  //!
  //! @return Iterator pointing to the inserted element.
  template<typename... Args>
  typename ListType::iterator emplace(const T1& order_item, Args&&... args) {
    auto order_it = order_items_.begin();
    auto insert_order_it = order_items_.before_begin();
    const auto order_end = order_items_.end();
    auto insert_item_it = items_.before_begin();

    while (order_it != order_end && comp_(*order_it, order_item)) {
      ++insert_item_it;
      ++insert_order_it;
      ++order_it;
    }
    order_items_.insert_after(insert_order_it, order_item);
    return items_.emplace_after(insert_item_it, std::forward<Args>(args)...);
  }

  //! Get the list of actual items.
  const std::forward_list<T2>& items() const noexcept {
    return items_;
  }

  //! Get the list of actual items.
  std::forward_list<T2>& items() noexcept {
    return items_;
  }

  typename ListType::iterator begin() noexcept {
    return items_.begin();
  }

  typename ListType::iterator end() noexcept {
    return items_.end();
  }

  typename ListType::const_iterator begin() const noexcept {
    return items_.begin();
  }

  typename ListType::const_iterator end() const noexcept {
    return items_.end();
  }

  typename ListType::const_iterator cbegin() const noexcept {
    return items_.cbegin();
  }

  typename ListType::const_iterator cend() const noexcept {
    return items_.cend();
  }

 private:
  Comparator comp_;
  std::forward_list<T1> order_items_;
  std::forward_list<T2> items_;
};
}  // namespace utility
}  // namespace pense

#endif  // CONTAINER_UTILITY_HPP_
