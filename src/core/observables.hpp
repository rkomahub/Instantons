#pragma once
#include "utils/io.hpp"
#include <vector>

/**
 * @file observables.hpp
 * @brief Measurement routines for Euclidean lattice configurations.
 *
 * Provides action, moments, and correlation functions
 * used in Monte Carlo and semiclassical analyses.
 */

/**
 * @brief Compute the Euclidean lattice action.
 *
 * Evaluates the discretized action with periodic
 * boundary conditions.
 *
 * @param x    Euclidean path.
 * @param a    Lattice spacing.
 * @param eta  Double well parameter.
 *
 * @return Total action.
 */
double compute_action(const std::vector<double> &x, double a, double eta);

/**
 * @brief Compute the field moment \( \langle x^p \rangle \).
 *
 * @param x      Euclidean path.
 * @param power  Integer power \( p \).
 *
 * @return Configuration average of \( x^p \).
 */
double compute_moment(const std::vector<double> &x, int power);

/**
 * @brief Compute the power correlator \( C_p(\tau) \).
 *
 * Evaluates:
 *
 * \( C_p(\tau) = \langle x^p(0)x^p(\tau) \rangle \)
 *
 * with periodic boundary conditions.
 *
 * @param x      Euclidean path.
 * @param power  Integer power \( p \).
 *
 * @return Correlator values for all \( \tau \).
 */
std::vector<double> compute_correlator_power(const std::vector<double> &x,
                                             int power);

/**
 * @brief Compute the Euclidean two-point correlator.
 *
 * Evaluates:
 *
 * \( C(\tau) = \langle x(0)x(\tau) \rangle \)
 *
 * @param x Euclidean path.
 *
 * @return Correlator values for all \( \tau \).
 */
std::vector<double> compute_correlator(const std::vector<double> &x);