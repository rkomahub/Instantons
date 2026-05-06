#include "analysis/fig10_rilm.hpp"

#include "core/instanton.hpp"
#include "core/observables.hpp"
#include "models/rilm.hpp"
#include "utils/io.hpp"
#include "utils/parameters.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

void run_rilm_analysis(std::mt19937 &gen) {
  const int trials = 50;
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;

  int n_inst = 4; // keep fixed for now
  if (n_inst % 2 != 0)
    ++n_inst; // safety

  std::vector<std::vector<double>> C1_data(trials), C2raw_data(trials),
      C2conn_data(trials), C3_data(trials);

  for (int t = 0; t < trials; ++t) {
    auto path = generate_rilm_path(N, a, eta, n_inst, gen);

    if (t == 0) {
      auto [mn, mx] = std::minmax_element(path.begin(), path.end());
      std::cout << "trial0 path min=" << *mn << " max=" << *mx << "\n";
      save_path_to_csv(path, "data/fig12_rilm_path.csv", a);
    }

    const double mean_x2 = compute_moment(path, 2);

    auto C1 = compute_correlator_power(path, 1);
    auto C2raw = compute_correlator_power(path, 2);
    auto C3 = compute_correlator_power(path, 3);

    auto C2conn = C2raw;
    for (int i = 0; i < N; ++i)
      C2conn[i] -= mean_x2 * mean_x2;

    // store (now all vectors have size N)
    C1_data[t] = std::move(C1);
    C2raw_data[t] = std::move(C2raw);
    C2conn_data[t] = std::move(C2conn);
    C3_data[t] = std::move(C3);
  }

  // ------------------------------
  // Export per-trial file (AFTER loop)
  // ------------------------------
  std::ofstream outT("data/fig10_trials_correlators.csv");
  outT << "trial,tau,C1,C2raw,C2conn,C3\n";
  for (int t = 0; t < trials; ++t) {
    for (int i = 0; i < N; ++i) {
      outT << t << "," << i * a << "," << C1_data[t][i] << ","
           << C2raw_data[t][i] << "," << C2conn_data[t][i] << ","
           << C3_data[t][i] << "\n";
    }
  }

  // ------------------------------
  // Compute ensemble averages (AFTER loop)
  // ------------------------------
  std::vector<double> C1_mean(N, 0), C1_err(N, 0);
  std::vector<double> C2raw_mean(N, 0), C2raw_err(N, 0);
  std::vector<double> C2conn_mean(N, 0), C2conn_err(N, 0);
  std::vector<double> C3_mean(N, 0), C3_err(N, 0);

  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < trials; ++t) {
      C1_mean[i] += C1_data[t][i];
      C2raw_mean[i] += C2raw_data[t][i];
      C2conn_mean[i] += C2conn_data[t][i];
      C3_mean[i] += C3_data[t][i];
    }
    C1_mean[i] /= trials;
    C2raw_mean[i] /= trials;
    C2conn_mean[i] /= trials;
    C3_mean[i] /= trials;

    double v1 = 0, v2r = 0, v2c = 0, v3 = 0;
    for (int t = 0; t < trials; ++t) {
      v1 += std::pow(C1_data[t][i] - C1_mean[i], 2);
      v2r += std::pow(C2raw_data[t][i] - C2raw_mean[i], 2);
      v2c += std::pow(C2conn_data[t][i] - C2conn_mean[i], 2);
      v3 += std::pow(C3_data[t][i] - C3_mean[i], 2);
    }
    if (trials > 1) {
      v1 /= (trials - 1);
      v2r /= (trials - 1);
      v2c /= (trials - 1);
      v3 /= (trials - 1);
    }
    C1_err[i] = std::sqrt(v1 / trials);
    C2raw_err[i] = std::sqrt(v2r / trials);
    C2conn_err[i] = std::sqrt(v2c / trials);
    C3_err[i] = std::sqrt(v3 / trials);
  }

  std::ofstream out("data/fig10_rilm_correlators.csv");
  out << "tau,C1,C1_err,C2raw,C2raw_err,C2conn,C2conn_err,C3,C3_err\n";
  for (int i = 0; i < N; ++i) {
    out << i * a << "," << C1_mean[i] << "," << C1_err[i] << ","
        << C2raw_mean[i] << "," << C2raw_err[i] << "," << C2conn_mean[i] << ","
        << C2conn_err[i] << "," << C3_mean[i] << "," << C3_err[i] << "\n";
  }

  std::cout << "[✓] Fig.10 RILM ensemble generated\n";
}