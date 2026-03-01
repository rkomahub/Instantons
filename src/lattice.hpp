#pragma once
#include <random>
#include <vector>

/**
 * @file lattice.hpp
 * @brief Representation of a one-dimensional Euclidean lattice configuration.
 *
 * The lattice stores the discretized path
 *
 *     x(τ_i),  i = 0, ..., N-1
 *
 * with periodic boundary conditions. It represents one configuration
 * sampled from the Euclidean path integral.
 */

class Lattice {
public:
  /**
   * @brief Construct a lattice configuration.
   *
   * @param N          Number of lattice points.
   * @param eta        Parameter defining the minima of the double well
   * potential.
   * @param hot_start  If true, initialize randomly in [-eta, eta].
   *                   If false, initialize in the classical vacuum x = +eta.
   * @param gen        Random number generator (injected).
   *
   * The hot start corresponds to a disordered configuration,
   * while the cold start corresponds to a classical minimum.
   */
  Lattice(int N, double eta, bool hot_start, std::mt19937 &gen);

  double &operator[](int i);
  const double &operator[](int i) const;
  int size() const;
  const std::vector<double> &get_path() const;

private:
  std::vector<double> x;
};