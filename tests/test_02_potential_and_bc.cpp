#include "parameters.hpp"
#include "potential.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

static double total_action(const std::vector<double> &x) {
  const int N = static_cast<int>(x.size());
  double S = 0.0;

  for (int i = 0; i < N; ++i) {
    const int ip = (i + 1) % N;
    S += std::pow(x[ip] - x[i], 2) / (4.0 * params::a);
    S += params::a * potential(x[i], params::eta);
  }

  return S;
}

int main() {

  // ---------- TEST 1: Potential symmetry V(x) = V(-x) ----------
  {
    const double xs[] = {-2.0, -1.0, -0.3, 0.0, 0.2, 1.1, 2.7};
    for (double x : xs) {
      const double v1 = potential(x, params::eta);
      const double v2 = potential(-x, params::eta);
      assert(std::abs(v1 - v2) < 1e-14);
    }
  }

  // ---------- TEST 2: Periodic boundary consistency ----------
  // The action must be invariant under a cyclic shift of the path.
  {
    const int N = 16;
    std::vector<double> x(N);
    for (int i = 0; i < N; ++i)
      x[i] = 0.1 * i - 0.3;

    std::vector<double> y(N);
    for (int i = 0; i < N; ++i)
      y[i] = x[(i + 1) % N]; // cyclic shift

    const double Sx = total_action(x);
    const double Sy = total_action(y);

    assert(std::abs(Sx - Sy) < 1e-12);
  }

  std::cout << "[✓] Potential + boundary-condition tests passed.\n";
  return 0;
}