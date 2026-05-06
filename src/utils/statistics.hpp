#pragma once

#include <vector>

/**
 * @file statistics.hpp
 * @brief Basic statistical utilities for Monte Carlo measurements.
 */

double mean(const std::vector<double> &values);

double sample_variance(const std::vector<double> &values, double mean_value);

double standard_error(const std::vector<double> &values);