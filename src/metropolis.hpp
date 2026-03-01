#pragma once

#include "lattice.hpp"
#include <random>

/**
 * @file metropolis.hpp
 * @brief Metropolis algorithm for Euclidean path integral sampling.
 *
 * This class performs Monte Carlo updates of a lattice configuration
 * according to the probability distribution
 *
 *     P[x] ∝ exp(−S_E[x])
 *
 * where S_E is the discretized Euclidean action of the double well model.
 *
 * The algorithm ensures detailed balance and ergodicity.
 */

class Metropolis {
public:
  /**
   * @brief Construct Metropolis evolver for a given lattice configuration.
   *
   * @param lattice Reference to configuration to be updated.
   *
   * The lattice is updated in-place.
   */
  Metropolis(Lattice &lattice);

  /**
   * @brief Perform one full Metropolis sweep.
   *
   * A sweep consists of attempting an update at every lattice site:
   *
   *     x_i → x_i + δx
   *
   * where δx is drawn from a Gaussian distribution.
   *
   * The proposal is accepted with probability
   *
   *     min(1, exp(−ΔS))
   *
   * ensuring sampling according to exp(−S_E).
   */
  void step();

  /// Perform n_sweeps of cooling (accept only action-decreasing updates).
  void cool(int n_sweeps);

private:
  /**
   * @brief Reference to lattice configuration.
   */
  Lattice &x;

  /**
   * @brief Random number generator used for updates.
   */
  std::mt19937 gen;
};