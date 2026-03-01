#pragma once
#include <vector>

/**
 * @file rilm.hpp
 * @brief Random Instanton Liquid Model (RILM).
 *
 * In the random instanton approximation, a configuration is constructed
 * as a superposition of independently distributed instantons and
 * anti-instantons with random positions.
 *
 * The multi-instanton ansatz used is:
 *
 *     x(τ) = η [ Σ_j Q_j tanh( ω/2 (τ − τ_j) ) − 1 ]
 *
 * where:
 *   Q_j = ±1 (instanton / anti-instanton)
 *   τ_j = random position in Euclidean time
 *   ω = 4 η
 */

/**
 * @brief Generate a random instanton liquid configuration.
 *
 * Instantons and anti-instantons are placed randomly
 * in Euclidean time without interactions.
 *
 * The total number of tunneling events is fixed to n_inst.
 *
 * @param N        Number of lattice sites.
 * @param a        Lattice spacing.
 * @param eta      Position of potential minima.
 * @param n_inst   Total number of instantons + anti-instantons.
 *
 * @return Discretized multi-instanton configuration.
 */
std::vector<double> generate_rilm_path(int N, double a, double eta, int n_inst);