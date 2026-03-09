#include "observables.hpp"
#include "potential.hpp"
#include <cmath>
#include <fstream>

double compute_action(const std::vector<double> &x, double a, double eta) {
  const int N = x.size();
  double S = 0.0;

  for (int i = 0; i < N; i++) {
    int im = (i - 1 + N) % N;

    double kinetic = std::pow(x[i] - x[im], 2) / (4.0 * a);

    double potential_term = a * potential(x[i], eta);

    S += kinetic + potential_term;
  }

  return S;
}

double compute_moment(const std::vector<double> &x, int power) {
  const int N = static_cast<int>(x.size());
  double sum = 0.0;

  for (int i = 0; i < N; ++i) {
    sum += std::pow(x[i], power);
  }

  return sum / static_cast<double>(N);
}

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
      sum += std::pow(x[i], power) * std::pow(x[j], power);
    }

    correlator[tau] = sum / static_cast<double>(N);
  }

  return correlator;
}

std::vector<double> compute_correlator(const std::vector<double> &x) {
  return compute_correlator_power(x, 1);
}

void save_correlator_to_csv(const std::vector<double> &correlator, double a,
                            const std::string &filename) {
  std::ofstream out(filename);

  out << "tau,C(tau)\n";

  for (size_t i = 0; i < correlator.size(); ++i) {
    out << i * a << "," << correlator[i] << "\n";
  }
}