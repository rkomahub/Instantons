#pragma once
#include <vector>

/**
 * @file iilm.hpp
 * @brief Interacting Instanton Liquid Model (IILM).
 *
 * In contrast to the Random Instanton Liquid Model (RILM),
 * the IILM includes short-range repulsive interactions
 * between instantons and anti-instantons.
 *
 * This prevents arbitrarily close instanton–anti-instanton pairs,
 * which would correspond to non-semiclassical configurations.
 *
 * The configuration is constructed using the same sum ansatz as RILM,
 * but with a minimal separation constraint.
 */

/**
 * @brief Generate an interacting instanton liquid configuration.
 *
 * Instantons and anti-instantons are placed randomly in Euclidean time,
 * subject to a minimal separation constraint:
 *
 *     |τ_i − τ_j| ≥ tau_core
 *
 * This mimics the short-range repulsive core in semiclassical
 * instanton ensembles.
 *
 * @param N         Number of lattice sites.
 * @param a         Lattice spacing.
 * @param eta       Position of potential minima.
 * @param n_inst    Target number of instantons + anti-instantons.
 * @param tau_core  Minimal allowed separation between instantons.
 *
 * @return Discretized interacting multi-instanton configuration.
 */
std::vector<double> generate_iilm_path(int N, double a, double eta, int n_inst,
                                       double tau_core);