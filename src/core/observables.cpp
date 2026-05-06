#include "core/observables.hpp"
#include "core/potential.hpp"
#include "utils/periodic.hpp"

#include <cmath>

namespace {

// Helper for powers of lattice field values.
double integer_power(double x, int power) { return std::pow(x, power); }

} // namespace

// Compute the discretized Euclidean action of one path.
double compute_action(const std::vector<double> &x, double a, double eta) {
  const int N = x.size();
  double S = 0.0;

  for (int i = 0; i < N; ++i) {
    const int im = periodic_prev(i, N);

    const double dx = x[i] - x[im];
    const double kinetic = dx * dx / (4.0 * a);
    const double potential_term = a * potential(x[i], eta);

    S += kinetic + potential_term;
  }

  return S;
}

// Compute the lattice average of x^power.
double compute_moment(const std::vector<double> &x, int power) {
  const int N = static_cast<int>(x.size());
  double sum = 0.0;

  for (int i = 0; i < N; ++i) {
    sum += integer_power(x[i], power);
  }

  return sum / static_cast<double>(N);
}

// Compute <x^power(0) x^power(tau)> on a periodic lattice.
std::vector<double> compute_correlator_power(const std::vector<double> &x,
                                             int power) {
  const int N = static_cast<int>(x.size());
  std::vector<double> correlator(N, 0.0);

  // For each Euclidean time separation τ,
  // average over all starting points i (translational invariance).
  for (int tau = 0; tau < N; ++tau) {
    double sum = 0.0;

    for (int i = 0; i < N; ++i) {
      const int j = (i + tau) % N; // periodic boundary conditions
      const double xi = integer_power(x[i], power);
      const double xj = integer_power(x[j], power);
      sum += xi * xj;
    }

    correlator[tau] = sum / static_cast<double>(N);
  }

  return correlator;
}

// Convenience wrapper for the coordinate correlator <x(0)x(tau)>.
std::vector<double> compute_correlator(const std::vector<double> &x) {
  return compute_correlator_power(x, 1);
}