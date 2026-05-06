#pragma once

/**
 * @file parameters.hpp
 * @brief Global parameters for the Euclidean lattice simulation.
 *
 * Defines lattice geometry, Monte Carlo settings,
 * and double well potential parameters.
 */

namespace params {

/**
 * @brief Number of Euclidean lattice sites.
 *
 * The Euclidean time extent is
 * \( \beta = Na \).
 */
constexpr int N = 800;

/**
 * @brief Euclidean lattice spacing.
 */
constexpr double a = 0.05;

/**
 * @brief Double well parameter.
 *
 * The potential minima are located at
 * \( x = \pm \eta \).
 */
inline double eta = 1.4;

/**
 * @brief Number of Metropolis sweeps.
 */
constexpr int sweeps = 100000;

/**
 * @brief Proposal width for Metropolis updates.
 */
constexpr double dx_width = 0.5;

/**
 * @brief Proposal width used during cooling.
 *
 * Cooling accepts only action-decreasing updates.
 */
constexpr double dx_width_cool = 0.05;

/**
 * @brief Initial configuration type.
 *
 * true  -> hot start
 * false -> cold start
 */
constexpr bool hot_start = true;

} // namespace params