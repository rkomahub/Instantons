#pragma once
#include <vector>

/**
 * @file heating.hpp
 * @brief Addition of Gaussian quantum fluctuations to classical configurations.
 *
 * Heating is the inverse of cooling:
 * starting from a classical (instanton) configuration,
 * short-distance Gaussian fluctuations are added in order
 * to mimic quantum fluctuations.
 *
 * This corresponds to sampling small fluctuations around
 * semiclassical backgrounds.
 */

/**
 * @brief Add Gaussian fluctuations to a path configuration.
 *
 * Each lattice site is modified as:
 *
 *     x_i → x_i + N(0, sigma)
 *
 * where N(0, sigma) is a Gaussian random variable.
 *
 * This approximates quantum fluctuations around a classical path.
 *
 * @param path   Configuration to be modified (in-place).
 * @param sigma  Standard deviation of Gaussian noise.
 */
void apply_gaussian_fluctuations(std::vector<double> &path, double sigma);