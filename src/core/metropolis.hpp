#pragma once

#include "core/lattice.hpp"
#include <random>

/**
 * @file metropolis.hpp
 * @brief Metropolis updates and cooling for Euclidean lattice paths.
 *
 * Samples configurations with probability
 * \( P[x] \propto e^{-S_E[x]} \)
 * for the double well Euclidean action.
 */

class Metropolis {
public:
  /**
   * @brief Construct a Metropolis evolver.
   *
   * The lattice configuration is updated in-place.
   *
   * @param lattice Lattice configuration.
   * @param gen     Random number generator.
   */
  Metropolis(Lattice &lattice, std::mt19937 &gen);

  /**
   * @brief Perform one full Metropolis sweep.
   */
  void step();

  /**
   * @brief Perform cooling sweeps.
   *
   * Cooling accepts only action-decreasing updates.
   *
   * @param n_sweeps Number of cooling sweeps.
   */
  void cool(int n_sweeps);

private:
  Lattice &x;

  /// Random number generator (non-owning reference).
  std::mt19937 &gen;
};