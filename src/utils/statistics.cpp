#include "utils/statistics.hpp"

#include <cmath>
#include <stdexcept>

// Compute the arithmetic mean of a data sample.
double mean(const std::vector<double> &values) {
  if (values.empty()) {
    throw std::runtime_error("mean: empty input vector");
  }

  double sum = 0.0;
  for (double x : values) {
    sum += x;
  }

  return sum / static_cast<double>(values.size());
}

// Compute the unbiased sample variance around a known mean.
double sample_variance(const std::vector<double> &values, double mean_value) {
  if (values.size() < 2) {
    return 0.0;
  }

  double var = 0.0;
  for (double x : values) {
    const double dx = x - mean_value;
    var += dx * dx;
  }

  return var / static_cast<double>(values.size() - 1);
}

// Compute the statistical standard error of the sample mean.
double standard_error(const std::vector<double> &values) {
  if (values.empty()) {
    throw std::runtime_error("standard_error: empty input vector");
  }

  const double m = mean(values);
  const double var = sample_variance(values, m);

  return std::sqrt(var / static_cast<double>(values.size()));
}