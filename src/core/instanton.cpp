#include "core/instanton.hpp"

// Count tunneling events through sign changes of the Euclidean path.
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

  // Each crossing corresponds to one instanton or anti-instanton event.
  return count;
}