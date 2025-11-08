#include "iilm.hpp"
#include <cmath>
#include <random>

std::vector<double> generate_iilm_path(int N, double a, double eta, int n_inst,
                                       double tau_core) {
  std::vector<double> x(N, 0.0);
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> dist_tau(0.0, N * a);
  std::uniform_int_distribution<int> dist_sign(0, 1);

  double omega = 4.0 * eta;

  std::vector<double> tau_i;
  std::vector<int> Q;

  int attempts = 0;
  while (tau_i.size() < static_cast<size_t>(n_inst) && attempts < 10 * n_inst) {
    double trial_tau = dist_tau(gen);
    bool too_close = false;

    for (double existing : tau_i) {
      if (std::abs(trial_tau - existing) < tau_core) {
        too_close = true;
        break;
      }
    }

    if (!too_close) {
      tau_i.push_back(trial_tau);
      Q.push_back(dist_sign(gen) == 0 ? -1 : 1);
    }
    ++attempts;
  }

  for (int i = 0; i < N; ++i) {
    double tau = i * a;
    double sum = 0.0;

    for (size_t j = 0; j < tau_i.size(); ++j) {
      sum += Q[j] * std::tanh(0.5 * omega * (tau - tau_i[j]));
    }

    x[i] = eta * (sum - 1.0);
  }

  return x;
}