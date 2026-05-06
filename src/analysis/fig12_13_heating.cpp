#include "analysis/fig12_13_heating.hpp"

#include "core/observables.hpp"
#include "models/heating.hpp"
#include "models/rilm.hpp"
#include "utils/io.hpp"
#include "utils/parameters.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

void run_heated_rilm_analysis(std::mt19937 &gen) {
  const int trials = 100; // before 50
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;

  // --- RILM parameters (match Fig.10 choice)
  int n_inst = 4;
  if (n_inst % 2 != 0)
    ++n_inst;

  // --- heating parameters (Fig.12/13)
  const int n_white_sweeps = 10;
  const double fluct_sigma = 0.2;

  const int n_heat_sweeps = 300;    // metro sweeps
  const double dx_width_heat = 0.3; // tune acceptance

  auto Vdd = [eta](double xc) { return 12.0 * xc * xc - 4.0 * eta * eta; };
  std::normal_distribution<double> dx_dist(0.0, dx_width_heat);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  // --- storage for Fig.13 correlators (heated ensemble)
  std::vector<std::vector<double>> C1_data(trials), C2raw_data(trials),
      C2conn_data(trials), C3_data(trials);

  long long n_prop_tot = 0, n_acc_tot = 0;

  for (int t = 0; t < trials; ++t) {

    // 1) fresh backbone per trial
    auto x_cl = generate_rilm_path(N, a, eta, n_inst, gen);

    // 2) for t==0 export Fig.12 typical paths (white + metro) from this
    // backbone
    if (t == 0) {
      save_path_to_csv(x_cl, "data/fig12_rilm_path.csv", a);

      // (A) white-noise heated typical path
      auto x_white = x_cl;
      for (int s = 0; s < n_white_sweeps; ++s) {
        apply_gaussian_fluctuations(x_white, fluct_sigma, gen);
      }
      save_path_to_csv(x_white, "data/fig12_heated_rilm_path.csv", a);

      // (B) metro heated typical path
      auto x_metro = x_cl;
      for (int sweep = 0; sweep < n_heat_sweeps; ++sweep) {
        for (int i = 0; i < N; ++i) {
          const int ip = (i + 1) % N;
          const int im = (i - 1 + N) % N;

          const double x_old = x_metro[i];
          const double kappa = Vdd(x_cl[i]);
          const double dx0 = x_old - x_cl[i];

          const double S_old = (std::pow(x_old - x_metro[im], 2) +
                                std::pow(x_metro[ip] - x_old, 2)) /
                                   (4.0 * a) +
                               0.5 * a * kappa * dx0 * dx0;

          const double x_new = x_old + dx_dist(gen);
          const double dx1 = x_new - x_cl[i];

          const double S_new = (std::pow(x_new - x_metro[im], 2) +
                                std::pow(x_metro[ip] - x_new, 2)) /
                                   (4.0 * a) +
                               0.5 * a * kappa * dx1 * dx1;

          const double dS = S_new - S_old;

          ++n_prop_tot;
          if (u01(gen) <= std::exp(-dS)) {
            x_metro[i] = x_new;
            ++n_acc_tot;
          }
        }
      }
      save_path_to_csv(x_metro, "data/fig12_heated_metro_rilm_path.csv", a);
    }

    // 3) metro-heated configuration for correlators (Fig.13 ensemble)
    std::vector<double> path = x_cl;
    for (int sweep = 0; sweep < n_heat_sweeps; ++sweep) {
      for (int i = 0; i < N; ++i) {
        const int ip = (i + 1) % N;
        const int im = (i - 1 + N) % N;

        const double x_old = path[i];
        const double kappa = Vdd(x_cl[i]);
        const double dx0 = x_old - x_cl[i];

        const double S_old =
            (std::pow(x_old - path[im], 2) + std::pow(path[ip] - x_old, 2)) /
                (4.0 * a) +
            0.5 * a * kappa * dx0 * dx0;

        const double x_new = x_old + dx_dist(gen);
        const double dx1 = x_new - x_cl[i];

        const double S_new =
            (std::pow(x_new - path[im], 2) + std::pow(path[ip] - x_new, 2)) /
                (4.0 * a) +
            0.5 * a * kappa * dx1 * dx1;

        const double dS = S_new - S_old;

        ++n_prop_tot;
        if (u01(gen) <= std::exp(-dS)) {
          path[i] = x_new;
          ++n_acc_tot;
        }
      }
    }

    // 4) correlators (same as Fig.10/4)
    const double mean_x2 = compute_moment(path, 2);

    auto C1 = compute_correlator_power(path, 1);
    auto C2raw = compute_correlator_power(path, 2);
    auto C3 = compute_correlator_power(path, 3);

    auto C2conn = C2raw;
    for (int i = 0; i < N; ++i)
      C2conn[i] -= mean_x2 * mean_x2;

    C1_data[t] = std::move(C1);
    C2raw_data[t] = std::move(C2raw);
    C2conn_data[t] = std::move(C2conn);
    C3_data[t] = std::move(C3);
  }

  // ---- export per-trial for jackknife (Fig.13)
  {
    std::ofstream outT("data/fig13_trials_correlators.csv");
    outT << "trial,tau,C1,C2raw,C2conn,C3\n";
    for (int t = 0; t < trials; ++t) {
      for (int i = 0; i < N; ++i) {
        outT << t << "," << i * a << "," << C1_data[t][i] << ","
             << C2raw_data[t][i] << "," << C2conn_data[t][i] << ","
             << C3_data[t][i] << "\n";
      }
    }
  }

  // ---- ensemble averages + stderr (Fig.13)
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

  {
    std::ofstream out("data/fig13_heated_rilm_correlators.csv");
    out << "tau,C1,C1_err,C2raw,C2raw_err,C2conn,C2conn_err,C3,C3_err\n";
    for (int i = 0; i < N; ++i) {
      out << i * a << "," << C1_mean[i] << "," << C1_err[i] << ","
          << C2raw_mean[i] << "," << C2raw_err[i] << "," << C2conn_mean[i]
          << "," << C2conn_err[i] << "," << C3_mean[i] << "," << C3_err[i]
          << "\n";
    }
  }

  const double acc_rate =
      (n_prop_tot > 0) ? double(n_acc_tot) / double(n_prop_tot) : 0.0;

  std::cout << "[✓] Fig.12 paths written (t=0) and Fig.13 heated-RILM "
               "correlators written (trials="
            << trials << ", metro acceptance=" << acc_rate << ")\n";
}