#include "instanton.hpp"

/**
 * @file instanton.cpp
 * @brief Identification of instantons via zero crossings.
 *
 * In one-dimensional quantum mechanics, instantons correspond to
 * tunneling events between the two classical vacua ±eta.
 *
 * On the Euclidean lattice, a tunneling event is identified
 * by a sign change of the field between neighboring sites:
 *
 *     x_i * x_{i+1} < 0
 *
 * Each such crossing corresponds to either:
 *   - an instanton      (− → +)
 *   - an anti-instanton (+ → −)
 *
 * Due to periodic boundary conditions, the lattice is treated
 * as a closed loop in Euclidean time.
 */

int count_zero_crossings(const std::vector<double> &path) {
  int count = 0;
  const int N = path.size();

  for (int i = 0; i < N; ++i) {
    const int j = (i + 1) % N; // periodic boundary condition

    // Detect sign change between neighboring lattice sites.
    if ((path[i] < 0 && path[j] > 0) || (path[i] > 0 && path[j] < 0)) {
      ++count;
    }
  }

  return count;
}