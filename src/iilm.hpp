#pragma once
#include <random>
#include <vector>

/**
 * @file iilm.hpp
 * @brief Interacting Instanton Liquid Model (IILM) utilities.
 *
 * For Fig.17 we need to evolve the instanton collective coordinates
 * (positions and charges) across many configurations.
 */

struct IILMConfig {
  std::vector<double> tau; // positions in [0,beta)
  std::vector<int> Q;      // +1 instanton, -1 anti-instanton
};

/**
 * @brief Generate an initial IILM collective-coordinate configuration.
 *
 * Places n_inst objects at random times in [0,beta), with periodic minimal
 * separation >= tau_core. Charges are assigned randomly (±1).
 *
 * @param n_inst     Total number of instantons + anti-instantons.
 * @param beta       Euclidean time extent.
 * @param tau_core   Minimal allowed periodic separation.
 * @param gen        RNG (injected).
 */
IILMConfig generate_iilm_config(int n_inst, double beta, double tau_core,
                                std::mt19937 &gen);

/**
 * @brief Build the sum-ansatz lattice path x_i from an IILMConfig.
 *
 * Uses the standard multi-(anti)instanton sum ansatz:
 *   x(tau) = eta * ( sum_j Q_j tanh(omega/2 (tau - tau_j)) - 1 )
 * with periodic nearest-image convention for (tau - tau_j).
 *
 * @param N     Number of lattice sites.
 * @param a     Lattice spacing.
 * @param eta   Double-well parameter.
 * @param cfg   Collective-coordinate configuration.
 */
std::vector<double> build_iilm_path(int N, double a, double eta,
                                    const IILMConfig &cfg);

/**
 * @brief Check core constraint: all pairs satisfy periodic_distance >=
 * tau_core.
 *
 * @param cfg       Configuration to test.
 * @param beta      Euclidean time extent.
 * @param tau_core  Minimal allowed periodic separation.
 */
bool iilm_core_ok(const IILMConfig &cfg, double beta, double tau_core);

/**
 * @brief One Metropolis proposal on the collective coordinates for Fig.17.
 *
 * Proposes to move one randomly chosen object:
 *   tau_j -> wrap(tau_j + d)
 * where d ~ Uniform(-step, step).
 *
 * This function only does the proposal+core-check; the driver will decide
 * acceptance based on the chosen effective action.
 *
 * @param cfg    Configuration to modify in-place if proposal passes core check.
 * @param beta   Euclidean time extent.
 * @param tau_core  Minimal allowed separation.
 * @param step   Proposal step size in time.
 * @param gen    RNG (injected).
 *
 * @return true if a *valid* (core-respecting) proposal was made and applied.
 *         false if rejected immediately by core constraint (cfg unchanged).
 */
bool propose_move_tau(IILMConfig &cfg, double beta, double tau_core,
                      double step, std::mt19937 &gen);