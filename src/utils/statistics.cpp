#include "utils/statistics.hpp"

#include <cmath>
#include <stdexcept>

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

double standard_error(const std::vector<double> &values) {
  if (values.empty()) {
    throw std::runtime_error("standard_error: empty input vector");
  }

  const double m = mean(values);
  const double var = sample_variance(values, m);

  return std::sqrt(var / static_cast<double>(values.size()));
}