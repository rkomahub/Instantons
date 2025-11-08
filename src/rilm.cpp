#include "rilm.hpp"
#include <cmath>
#include <random>

std::vector<double> generate_rilm_path(int N, double a, double eta,
                                       int n_inst) {
  std::vector<double> x(N, 0.0);
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> dist_tau(0.0, N * a);
  std::uniform_int_distribution<int> dist_sign(0, 1);

  double omega = 4.0 * eta;

  std::vector<double> tau_i(n_inst);
  std::vector<int> Q(n_inst);

  for (int i = 0; i < n_inst; ++i) {
    tau_i[i] = dist_tau(gen);
    Q[i] = dist_sign(gen) == 0 ? -1 : 1;
  }

  for (int i = 0; i < N; ++i) {
    double tau = i * a;
    double sum = 0.0;

    for (int j = 0; j < n_inst; ++j) {
      sum += Q[j] * std::tanh(0.5 * omega * (tau - tau_i[j]));
    }

    x[i] = eta * (sum - 1.0);
  }

  return x;
}
