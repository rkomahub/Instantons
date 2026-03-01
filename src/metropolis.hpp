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
   * @param gen     Random number generator (injected).
   *
   * The lattice is updated in-place.
   */
  Metropolis(Lattice &lattice, std::mt19937 &gen);

  /**
   * @brief Perform one full Metropolis sweep.
   */
  void step();

  /// Perform n_sweeps of cooling (accept only action-decreasing updates).
  void cool(int n_sweeps);

private:
  Lattice &x;

  /**
   * @brief Random number generator used for updates.
   *
   * This generator is not owned by the class.
   */
  std::mt19937 &gen;
};