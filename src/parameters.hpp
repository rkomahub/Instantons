#pragma once

/**
 * @file parameters.hpp
 * @brief Global simulation parameters for the Euclidean lattice instanton
 * simulation.
 *
 * The system is the quantum mechanical double well potential
 *
 *     V(x) = (x^2 - eta^2)^2
 *
 * discretized on a one-dimensional Euclidean lattice of size N
 * with lattice spacing a and periodic boundary conditions.
 */

namespace params {

/**
 * @brief Number of lattice points in Euclidean time.
 *
 * The total Euclidean time extent is
 *
 *     beta = N * a
 *
 * In order to resolve instanton physics, one requires
 *     a << tau_osc
 *     beta >> tau_tunnel
 */
constexpr int N = 800;

/**
 * @brief Lattice spacing in Euclidean time.
 *
 * Must be small compared to the oscillation time scale:
 *     tau_osc ~ (4 eta)^{-1}
 */
constexpr double a = 0.05;

/**
 * @brief Position of classical minima of the double well potential.
 *
 * The potential is
 *
 *     V(x) = (x^2 - eta^2)^2
 *
 * The classical instanton action is
 *
 *     S0 = 4 eta^3 / 3
 */
constexpr double eta = 1.4;

/**
 * @brief Number of full Metropolis sweeps.
 *
 * One sweep updates all lattice sites once.
 */
constexpr int sweeps = 100000;

/**
 * @brief Width of Gaussian proposal distribution.
 *
 * Trial updates are generated as
 *
 *     x_i → x_i + N(0, dx_width)
 *
 * Should be tuned to obtain ~50% acceptance rate.
 */
constexpr double dx_width = 0.5;

/**
 * @brief Width of Gaussian proposal distribution used during cooling.
 *
 * During cooling, updates are accepted only if they decrease the
 * Euclidean action. Unlike the Monte Carlo width (dx_width),
 * this parameter controls how aggressively the configuration is
 * smoothed toward a nearby classical solution.
 *
 * Smaller values make cooling behave closer to gradient descent,
 * efficiently suppressing UV fluctuations.
 *
 * Larger values make cooling more stochastic and may leave residual
 * oscillations even after many sweeps.
 *
 * Typical tuning range:
 *     0.05 – 0.20
 *
 * This parameter does NOT affect full Monte Carlo sampling.
 */
constexpr double dx_width_cool = 0.05;

/**
 * @brief Initial configuration selector.
 *
 * true  → hot start (random in [-eta, eta])
 * false → cold start (x_i = +eta)
 */
constexpr bool hot_start = true;

} // namespace params