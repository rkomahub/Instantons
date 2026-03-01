#include "observables.hpp"
#include <fstream>

/**
 * @file observables.cpp
 * @brief Implementation of Euclidean two-point correlators.
 *
 * Given a configuration x_i ≡ x(τ_i) sampled from the Euclidean path integral,
 * this module computes the translationally averaged correlator
 *
 *     C(τ) = (1/N) Σ_i x_i x_{i+τ}
 *
 * with periodic boundary conditions.
 *
 * In the large-τ limit:
 *
 *     C(τ) ~ exp(−(E1 − E0) τ)
 *
 * which allows extraction of the energy gap.
 */

std::vector<double> compute_correlator(const std::vector<double> &x) {
  const int N = x.size();
  std::vector<double> correlator(N, 0.0);

  // For each Euclidean time separation τ,
  // average over all starting points i (translational invariance).
  for (int tau = 0; tau < N; ++tau) {
    double sum = 0.0;

    for (int i = 0; i < N; ++i) {
      const int j = (i + tau) % N; // periodic boundary conditions
      sum += x[i] * x[j];
    }

    correlator[tau] = sum / static_cast<double>(N);
  }

  return correlator;
}

/**
 * @brief Save correlator to CSV file.
 *
 * The output format is:
 *
 *     tau, C(tau)
 *
 * where tau = i * a is the physical Euclidean time separation.
 */
void save_correlator_to_csv(const std::vector<double> &correlator, double a,
                            const std::string &filename) {
  std::ofstream out(filename);

  out << "tau,C(tau)\n";

  for (size_t i = 0; i < correlator.size(); ++i) {
    out << i * a << "," << correlator[i] << "\n";
  }
}