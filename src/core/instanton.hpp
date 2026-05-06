#pragma once
#include <vector>

/**
 * @file instanton.hpp
 * @brief Instanton identification utilities.
 *
 * In one-dimensional quantum mechanics, instantons correspond
 * to tunneling events between the two classical minima ±eta.
 *
 * On the Euclidean lattice, tunneling events can be identified
 * by zero-crossings of the path x(τ).
 */

/**
 * @brief Count zero crossings of a Euclidean path.
 *
 * A zero crossing occurs when consecutive lattice points satisfy
 *
 *     x_i * x_{i+1} < 0
 *
 * Each zero crossing corresponds to either:
 *   - an instanton  (− → + transition)
 *   - an anti-instanton (+ → − transition)
 *
 * The total number of instantons plus anti-instantons is
 *
 *     N_{I+A} = number of zero crossings.
 *
 * The instanton density is then
 *
 *     n_{I+A} = N_{I+A} / beta
 *
 * where beta = N * a is the Euclidean time extent.
 *
 * @param path  Discretized Euclidean path.
 * @return      Total number of zero crossings.
 */
int count_zero_crossings(const std::vector<double> &path);