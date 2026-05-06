#pragma once
#include <random>
#include <vector>

/**
 * @file heating.hpp
 * @brief Gaussian fluctuations around semiclassical configurations.
 *
 * Heating adds short-distance noise to classical instanton paths,
 * mimicking quantum fluctuations.
 */

/**
 * @brief Add Gaussian fluctuations to a path.
 *
 * Each lattice site is updated by Gaussian noise
 * with standard deviation \( \sigma \).
 *
 * @param path   Configuration modified in-place.
 * @param sigma  Width of the Gaussian fluctuations.
 * @param gen    Random number generator.
 */
void apply_gaussian_fluctuations(std::vector<double> &path, double sigma,
                                 std::mt19937 &gen);