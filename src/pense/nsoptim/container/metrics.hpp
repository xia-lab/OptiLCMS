//
//  metrics.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_CONTAINER_METRICS_HPP_
#define NSOPTIM_CONTAINER_METRICS_HPP_

#include <memory>
#include <string>
#include <forward_list>

#include "../config.hpp"

namespace nsoptim {
//! A named metric object.
template<typename T>
struct Metric {
  //! Construct a named metric object.
  Metric(const std::string& _name, const T _value) : name(_name), value(_value) {}

  //! Default move constructor.
  Metric(Metric&&) = default;
  //! Default move assignment.
  Metric& operator=(Metric&&) = default;
  //! Default copy constructor.
  Metric(const Metric&) = default;
  //! Default copy assignment.
  Metric& operator=(const Metric&) = default;

  std::string name;  //!< Name of the metric.
  T value;  //!< Value of the metric.
};

//! A list of Metrics
template <typename T> using MetricList = std::forward_list<Metric<T>>;

//! Structure to hold all supported metrics
struct MetricsContainer {
  MetricList<int> int_;
  MetricList<double> dbl_;
  MetricList<std::string> str_;
};

namespace _metrics_internal {
template<typename T>
class MetricsIterator : public std::iterator<std::forward_iterator_tag, const Metric<T>> {
  using list_iterator = typename MetricList<T>::const_iterator;

 public:
  MetricsIterator() noexcept {}
  // Create an iterator pointing to the start of the chained lists.
  MetricsIterator(const MetricList<T>& first, const MetricList<T>& second) noexcept
    : end_first_(first.cend()), begin_second_(second.cbegin()), hit_first_end_(first.cbegin() == end_first_),
      it_(hit_first_end_ ? begin_second_ : first.cbegin()) {}

  // Create an iterator pointing to the very end of the chained lists.
  MetricsIterator(const MetricList<T>& first, const MetricList<T>& second, std::true_type) noexcept
    : end_first_(first.cend()), begin_second_(second.cbegin()), hit_first_end_(true), it_(second.cend()) {}

  MetricsIterator(const MetricsIterator&) = default;
  MetricsIterator& operator=(const MetricsIterator&) = default;
  MetricsIterator(MetricsIterator&&) = default;
  MetricsIterator& operator=(MetricsIterator&&) = default;

  friend bool operator==(const MetricsIterator& lhs, const MetricsIterator& rhs) noexcept {
    return (lhs.hit_first_end_ == rhs.hit_first_end_) && (lhs.it_ == rhs.it_);
  }

  friend bool operator!=(const MetricsIterator& lhs, const MetricsIterator& rhs) noexcept {
    return (lhs.hit_first_end_ != rhs.hit_first_end_) || (lhs.it_ != rhs.it_);
  }

  MetricsIterator& operator++() noexcept {
    ++it_;
    if (it_ == end_first_ && !hit_first_end_) {
      it_ = begin_second_;
      hit_first_end_ = true;
    }
    return *this;
  }

  MetricsIterator operator++(int) noexcept {
    MetricsIterator other(*this);
    ++(*this);
    return other;
  }

  const Metric<T>& operator*() const noexcept {
    return *it_;
  }

  Metric<T> const * operator->() const noexcept {
    return &(*it_);
  }

 private:
  list_iterator end_first_;
  list_iterator begin_second_;
  bool hit_first_end_ = false;
  list_iterator it_;
};

//! Enumeration of states for metrics.
enum class HasMetrics {
  kNo,
  kMaybe,
  kYes
};

//! A collection of metrics
template<int>
class Metrics {
  //! A dummy iterator proxy over an empty set.
  template<typename T>
  class IteratorProxy {
   public:
    constexpr T const * begin() const noexcept { return nullptr; }
    constexpr T const * end() const noexcept { return nullptr; }
  };

 public:
  //! Construct a named collection of metrics.
  explicit Metrics(const std::string&) noexcept {}

  //! Get the name(space) of this Metrics object.
  //!
  //! @return The name(space).
  constexpr char const * name() const noexcept {
    return "";
  }

  //! Create a sub-collection of metrics and add it to this collection.
  //!
  //! @param name Name of the new sub-collection of metrics.
  //! @return Reference to the newly created metrics.
  Metrics& CreateSubMetrics(const std::string&) noexcept {
    return *this;
  }

  //! Add a sub-collection of Metrics to the collection.
  //!
  //! @param metrics A collection of metrics to be added to this metrics collection.
  void AddSubMetrics(const Metrics&) const noexcept {}

  //! Add a floating-point metric to the collection.
  //!
  //! @param name Name of the metric.
  //! @param value Value for the metric.
  void AddMetric(const std::string&, const double) const noexcept {}

  //! Add an integer metric to the collection.
  //!
  //! @param name Name of the metric.
  //! @param value Value for the metric.
  void AddMetric(const std::string&, const int) const noexcept {}

  //! Add a string metric to the collection.
  //!
  //! @param name Name of the metric.
  //! @param value Value for the metric.
  void AddMetric(const std::string&, const std::string&) const noexcept {}

  //! Add a floating-point detail to the collection.
  //!
  //! @param name Name of the detailed metric.
  //! @param value Value for the detailed metric.
  void AddDetail(const std::string&, const double) const noexcept {}

  //! Add an integer detail to the collection.
  //!
  //! @param name Name of the detailed metric.
  //! @param value Value for the detailed metric.
  void AddDetail(const std::string&, const int) const noexcept {}

  //! Add a string detail to the collection.
  //!
  //! @param name Name of the detailed metric.
  //! @param value Value for the detailed metric.
  void AddDetail(const std::string&, const std::string&) const noexcept {}

  //! Iterate over all floating-point metrics and details in this collection.
  //!
  //! @return A proxy providing a forward-iterator over @ref Metric "Metric<double>" objects.
  constexpr IteratorProxy<Metric<double>> DoubleMetrics() const noexcept {
    return IteratorProxy<Metric<double>>();
  }

  //! Iterate over all integer metrics and details in this collection.
  //!
  //! @return A proxy providing a forward-iterator over @ref Metric "Metric<int>" objects.
  constexpr IteratorProxy<Metric<int>> IntegerMetrics() const noexcept {
    return IteratorProxy<Metric<int>>();
  }

  //! Iterate over all string metrics and details in this collection.
  //!
  //! @return A proxy providing a forward-iterator over @ref Metric "Metric<std::string>" objects.
  constexpr IteratorProxy<Metric<std::string>> StringMetrics() const noexcept {
    return IteratorProxy<Metric<std::string>>();
  }

  //! Iterate over all nested metrics in this collection.
  //!
  //! @return A proxy providing a forward-iterator over Metrics objects.
  constexpr IteratorProxy<Metrics> SubMetrics() const noexcept {
    return IteratorProxy<Metrics>();
  }
};

template<>
class Metrics<1> {
 public:
  explicit Metrics(const std::string& name) noexcept : name_(name) {}

  Metrics(const Metrics& other) noexcept : name_(other.name_), metrics_(other.metrics_),
                                           has_metrics_(other.has_metrics_) {
    // We are not sure if there's actual data in the sub-metrics. Thin the empty ones out.
    for (auto&& sub_metric : other.sub_metrics_) {
      if (sub_metric.CheckMetrics() == HasMetrics::kYes) {
        sub_metrics_.push_front(sub_metric);
        has_metrics_ = HasMetrics::kYes;
        sub_metrics_.front().has_metrics_ = HasMetrics::kYes;
      }
    }
  }

  //! Get the name(space) of this Metrics object.
  //!
  //! @return the name(space).
  std::string name() const noexcept {
    return name_;
  }

  //! Create a sub-collection of metrics and add it to this collection.
  //!
  //! @param name the name of the new sub-collection of metrics.
  //! @return a reference to the newly created metrics.
  Metrics& CreateSubMetrics(const std::string& name) noexcept {
    has_metrics_ = HasMetrics::kMaybe;
    sub_metrics_.emplace_front(name);
    return sub_metrics_.front();
  }

  //! Add a sub-collection of metrics to the collection.
  //!
  //! @param metrics a collection of metrics to be added to this metrics collection.
  void AddSubMetrics(const Metrics& metrics) noexcept {
    auto sub_has_metrics = metrics.CheckMetrics();
    if (sub_has_metrics != HasMetrics::kNo) {
      has_metrics_ = sub_has_metrics;
      sub_metrics_.push_front(metrics);
      sub_metrics_.front().has_metrics_ = sub_has_metrics;
    }
  }

  //! Add a sub-collection of metrics to the collection.
  //!
  //! @param metrics a collection of metrics to be added to this metrics collection.
  void AddSubMetrics(Metrics&& metrics) noexcept {
    auto sub_has_metrics = metrics.CheckMetrics();
    if (sub_has_metrics != HasMetrics::kNo) {
      has_metrics_ = sub_has_metrics;
      sub_metrics_.emplace_front(std::move(metrics));
      sub_metrics_.front().has_metrics_ = sub_has_metrics;
    }
  }

  //! Add a metric to the collection.
  //!
  //! @param name the name of the metric.
  //! @param value the value for the metric.
  void AddMetric(const std::string& name, const double value) noexcept {
    has_metrics_ = HasMetrics::kYes;
    metrics_.dbl_.emplace_front(name, value);
  }

  //! Add a metric to the collection.
  //!
  //! @param name the name of the metric.
  //! @param value the value for the metric.
  void AddMetric(const std::string& name, const int value) noexcept {
    has_metrics_ = HasMetrics::kYes;
    metrics_.int_.emplace_front(name, value);
  }

  //! Add a metric to the collection.
  //!
  //! @param name the name of the metric.
  //! @param value the value for the metric.
  void AddMetric(const std::string& name, const std::string& value) noexcept {
    has_metrics_ = HasMetrics::kYes;
    metrics_.str_.emplace_front(name, value);
  }

  //! Add a detailed metric to the collection.
  //!
  //! @param name the name of the detailed metric.
  //! @param value the value for the detailed metric.
  void AddDetail(const std::string&, const double) const noexcept {}
  void AddDetail(const std::string&, const int) const noexcept {}
  void AddDetail(const std::string&, const std::string&) const noexcept {}

  const MetricList<double>& DoubleMetrics() const noexcept {
    return metrics_.dbl_;
  }
  const MetricList<int>& IntegerMetrics() const noexcept {
    return metrics_.int_;
  }
  const MetricList<std::string>& StringMetrics() const noexcept {
    return metrics_.str_;
  }

  const std::forward_list<Metrics>& SubMetrics() const noexcept {
    return sub_metrics_;
  }

 private:
  HasMetrics CheckMetrics() const noexcept {
    if (has_metrics_ != HasMetrics::kMaybe) {
      return has_metrics_;
    }

    // If we are not sure, check if the metrics are all empty.
    if (!metrics_.dbl_.empty() || !metrics_.int_.empty() || !metrics_.str_.empty()) {
      return HasMetrics::kYes;
    } else {
      // The metrics are all empty, but maybe one of the sub-metrics in non-empty.
      HasMetrics updated_has_metrics = HasMetrics::kNo;
      for (auto&& sub_metric : sub_metrics_) {
        auto sub_has_metrics = sub_metric.CheckMetrics();
        if (updated_has_metrics == HasMetrics::kYes) {
          // Stop after the first sub-metric has metrics.
          return HasMetrics::kYes;
        } else if (sub_has_metrics == HasMetrics::kMaybe) {
          updated_has_metrics = HasMetrics::kMaybe;
        }
      }
      return updated_has_metrics;
    }
  };

  std::string name_;
  MetricsContainer metrics_;
  std::forward_list<Metrics> sub_metrics_;
  HasMetrics has_metrics_ = HasMetrics::kNo;
};

template<>
class Metrics<2> {
 public:
  //! Proxy class to iterate over all the metrics of the same type.
  template<typename T>
  class MetricsIteratorProxy {
   public:
    using iterator = _metrics_internal::MetricsIterator<T>;

    //! Get an iterator pointing to the first element of the list.
    iterator begin() const noexcept { return iterator(metrics_, details_); }
    //! Get an iterator pointing past the last element of the list.
    iterator end() const noexcept { return iterator(metrics_, details_, std::true_type{}); }

   private:
    friend class Metrics;

    //! Constructor is only available to the Metrics<true> class.
    explicit MetricsIteratorProxy(const MetricList<T>& metrics, const MetricList<T>& details) noexcept
      : metrics_(metrics), details_(details) {}

    const MetricList<T>& metrics_;
    const MetricList<T>& details_;
  };

  explicit Metrics(const std::string& name) noexcept : name_(name) {}

  //! Get the name(space) of this Metrics object.
  //!
  //! @return the name(space).
  std::string name() const noexcept {
    return name_;
  }

  //! Create a sub-collection of metrics and add it to this collection.
  //!
  //! @param name the name of the new sub-collection of metrics.
  //! @return a reference to the newly created metrics.
  Metrics& CreateSubMetrics(const std::string& name) noexcept {
    sub_metrics_.emplace_front(name);
    return sub_metrics_.front();
  }

  //! Add a sub-collection of metrics to the collection.
  //!
  //! @param metrics a collection of metrics to be added to this metrics collection.
  void AddSubMetrics(const Metrics& metrics) noexcept {
    sub_metrics_.push_front(metrics);
  }

  //! Add a sub-collection of metrics to the collection.
  //!
  //! @param metrics a collection of metrics to be added to this metrics collection.
  void AddSubMetrics(Metrics&& metrics) noexcept {
    sub_metrics_.emplace_front(std::move(metrics));
  }

  //! Add a metric to the collection.
  //!
  //! @param name the name of the metric.
  //! @param value the value for the metric.
  void AddMetric(const std::string& name, const double value) noexcept {
    metrics_.dbl_.emplace_front(name, value);
  }

  //! Add a metric to the collection.
  //!
  //! @param name the name of the metric.
  //! @param value the value for the metric.
  void AddMetric(const std::string& name, const int value) noexcept {
    metrics_.int_.emplace_front(name, value);
  }

  //! Add a metric to the collection.
  //!
  //! @param name the name of the metric.
  //! @param value the value for the metric.
  void AddMetric(const std::string& name, const std::string& value) noexcept {
    metrics_.str_.emplace_front(name, value);
  }

  //! Add a detailed metric to the collection.
  //!
  //! @param name the name of the detailed metric.
  //! @param value the value for the detailed metric.
  void AddDetail(const std::string& name, const double value) noexcept {
    details_.dbl_.emplace_front(name, value);
  }

  //! Add a detailed metric to the collection.
  //!
  //! @param name the name of the detailed metric.
  //! @param value the value for the detailed metric.
  void AddDetail(const std::string& name, const int value) noexcept {
    details_.int_.emplace_front(name, value);
  }

  //! Add a detailed metric to the collection.
  //!
  //! @param name the name of the detailed metric.
  //! @param value the value for the detailed metric.
  void AddDetail(const std::string& name, const std::string& value) noexcept {
    details_.str_.emplace_front(name, value);
  }

  MetricsIteratorProxy<double> DoubleMetrics() const noexcept {
    return MetricsIteratorProxy<double>(metrics_.dbl_, details_.dbl_);
  }

  MetricsIteratorProxy<int> IntegerMetrics() const noexcept {
    return MetricsIteratorProxy<int>(metrics_.int_, details_.int_);
  }

  MetricsIteratorProxy<std::string> StringMetrics() const noexcept {
    return MetricsIteratorProxy<std::string>(metrics_.str_, details_.str_);
  }

  const std::forward_list<Metrics>& SubMetrics() const noexcept {
    return sub_metrics_;
  }

 private:
  std::string name_;
  MetricsContainer metrics_;
  MetricsContainer details_;
  std::forward_list<Metrics> sub_metrics_;
};
}  // namespace _metrics_internal
}  // namespace nsoptim

#endif  // NSOPTIM_CONTAINER_METRICS_HPP_
