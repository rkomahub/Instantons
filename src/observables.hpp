#pragma once
#include <string>
#include <vector>

/**
 * @file observables.hpp
 * @brief Measurement of physical observables from Euclidean configurations.
 *
 * This module provides tools to compute correlation functions
 * from a discretized Euclidean path x(τ_i).
 */

/**
 * @brief Compute the p-th moment of the field.
 *
 * Given a Euclidean configuration x_i ≡ x(τ_i), this function computes
 *
 *     ⟨ x^p ⟩ = (1/N) Σ_i x_i^p
 *
 * This quantity is needed for constructing connected correlators.
 *
 * For example, in the case p = 2:
 *
 *     C_2^conn(τ) = ⟨ x^2(0)x^2(τ) ⟩ − ⟨ x^2 ⟩^2
 *
 * @param x      Discretized Euclidean path.
 * @param power  Integer power p.
 * @return       The ensemble average ⟨ x^p ⟩ for this configuration.
 */
double compute_moment(const std::vector<double> &x, int power);

/**
 * @brief Compute the p-th power Euclidean two-point correlator.
 *
 * Given a configuration x_i ≡ x(τ_i), this function computes
 *
 *     C_p(τ) = (1/N) Σ_i x_i^p x_{i+τ}^p
 *
 * with periodic boundary conditions.
 *
 * Special cases:
 *
 *     p = 1  →  ⟨ x(0)x(τ) ⟩
 *     p = 2  →  ⟨ x^2(0)x^2(τ) ⟩
 *     p = 3  →  ⟨ x^3(0)x^3(τ) ⟩
 *
 * For p = 2, the connected correlator used in Fig. 4 is
 *
 *     C_2^conn(τ) = C_2(τ) − ⟨ x^2 ⟩^2
 *
 * @param x      Discretized Euclidean path.
 * @param power  Integer power p.
 * @return       Vector C_p(τ) for τ = 0, ..., N−1.
 */
std::vector<double> compute_correlator_power(const std::vector<double> &x,
                                             int power);

/**
 * @brief Compute the Euclidean two-point correlator.
 *
 * Given a configuration x_i ≡ x(τ_i), this function computes
 *
 *     C(τ) = (1/N) Σ_i x_i x_{i+τ}
 *
 * with periodic boundary conditions.
 *
 * The correlator is used to extract the energy gap via
 *
 *     C(τ) ~ exp(−(E1 − E0) τ)   as τ → ∞.
 *
 * @param x  Discretized Euclidean path.
 * @return   Vector C(τ) for τ = 0, ..., N−1.
 */
std::vector<double> compute_correlator(const std::vector<double> &x);

/**
 * @brief Save correlator data to CSV file.
 *
 * The output format is:
 *
 *     tau, C(tau)
 *
 * where tau = i * a.
 *
 * @param correlator  Correlator values.
 * @param a           Lattice spacing.
 * @param filename    Output file path.
 */
void save_correlator_to_csv(const std::vector<double> &correlator, double a,
                            const std::string &filename);