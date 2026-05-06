#pragma once
#include <random>
#include <vector>

/**
 * @file iilm.hpp
 * @brief Interacting Instanton Liquid Model utilities.
 *
 * Provides collective-coordinate configurations and
 * Monte Carlo updates for interacting instanton ensembles.
 */

struct IILMConfig {
  std::vector<double> tau; // instanton positions
  std::vector<int> Q;      // charges: +1 or -1
};

/**
 * @brief Generate an initial IILM configuration.
 *
 * Instanton positions are sampled in Euclidean time
 * with a minimal separation constraint.
 *
 * @param n_inst    Number of instantons and anti-instantons.
 * @param beta      Euclidean time extent.
 * @param tau_core  Minimal allowed separation.
 * @param gen       Random number generator.
 */
IILMConfig generate_iilm_config(int n_inst, double beta, double tau_core,
                                std::mt19937 &gen);

/**
 * @brief Build a lattice path from collective coordinates.
 *
 * Uses a multi-instanton sum ansatz with periodic boundary conditions.
 *
 * @param N     Number of lattice sites.
 * @param a     Lattice spacing.
 * @param eta   Double well parameter.
 * @param cfg   Instanton configuration.
 *
 * @return Euclidean lattice path.
 */
std::vector<double> build_iilm_path(int N, double a, double eta,
                                    const IILMConfig &cfg);

/**
 * @brief Check the instanton core constraint.
 *
 * @param cfg       Instanton configuration.
 * @param beta      Euclidean time extent.
 * @param tau_core  Minimal allowed separation.
 *
 * @return True if all separations satisfy the constraint.
 */
bool iilm_core_ok(const IILMConfig &cfg, double beta, double tau_core);

/**
 * @brief Propose a Monte Carlo move in collective-coordinate space.
 *
 * Moves one instanton position while enforcing
 * the core separation constraint.
 *
 * @param cfg       Configuration modified in-place.
 * @param beta      Euclidean time extent.
 * @param tau_core  Minimal allowed separation.
 * @param step      Proposal step size.
 * @param gen       Random number generator.
 *
 * @return True if the proposal satisfies the core constraint.
 */
bool propose_move_tau(IILMConfig &cfg, double beta, double tau_core,
                      double step, std::mt19937 &gen);