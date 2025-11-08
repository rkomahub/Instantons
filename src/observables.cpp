#include "observables.hpp"
#include <fstream>

std::vector<double> compute_correlator(const std::vector<double> &x) {
  int N = x.size();
  std::vector<double> correlator(N, 0.0);

  for (int tau = 0; tau < N; ++tau) {
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
      int j = (i + tau) % N;
      sum += x[i] * x[j];
    }
    correlator[tau] = sum / N;
  }
  return correlator;
}

void save_correlator_to_csv(const std::vector<double> &correlator, double a,
                            const std::string &filename) {
  std::ofstream out(filename);
  out << "tau,C(tau)\n";
  for (size_t i = 0; i < correlator.size(); ++i) {
    out << i * a << "," << correlator[i] << "\n";
  }
}
