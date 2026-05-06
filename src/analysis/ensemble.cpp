#include "analysis/ensemble.hpp"
#include "core/instanton.hpp"
#include "core/lattice.hpp"
#include "core/metropolis.hpp"
#include "core/observables.hpp"
#include "utils/parameters.hpp"
#include "utils/statistics.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// Generate many configurations and average observables over the ensemble.
void run_ensemble_average(int trials, bool cooled,
                          const std::string &output_prefix, std::mt19937 &gen) {

  const int N = params::N;
  const double a = params::a;
  const double beta = N * a;

  // For Fig. 3 (only uncooled)
  std::vector<double> all_positions;
  if (!cooled) {
    all_positions.reserve(static_cast<size_t>(trials) * static_cast<size_t>(N));
  }

  // Correlators per trial
  std::vector<std::vector<double>> C1_data(trials);
  std::vector<std::vector<double>> C2raw_data(trials);
  std::vector<std::vector<double>> C2conn_data(trials);
  std::vector<std::vector<double>> C3_data(trials);

  // Store one instanton count per configuration.
  std::vector<int> instanton_counts(trials);

  // ===============================
  // Generate ensemble
  // ===============================
  for (int t = 0; t < trials; ++t) {

    Lattice lattice(N, params::eta, /*hot_start=*/true, gen);
    Metropolis evolver(lattice, gen);

    // Generate one quantum configuration.
    for (int sweep = 0; sweep < params::sweeps; ++sweep) {
      evolver.step();
    }

    // Optionally remove short-distance fluctuations by cooling.
    if (cooled) {
      Lattice cooled_lat = lattice;
      Metropolis cooled_evolver(cooled_lat, gen);
      cooled_evolver.cool(200);
      lattice = cooled_lat;
    }

    const auto &path = lattice.get_path();

    // --- Fig. 3 positions ---
    if (!cooled) {
      all_positions.insert(all_positions.end(), path.begin(), path.end());
    }

    // --- Correlators ---
    auto C1 = compute_correlator_power(path, 1);
    auto C2raw = compute_correlator_power(path, 2);
    auto C3 = compute_correlator_power(path, 3);

    // Connected subtraction for x^2: C2conn = C2raw - <x^2>^2
    const double mean_x2 = compute_moment(path, 2);
    auto C2conn = C2raw;
    for (int i = 0; i < N; ++i) {
      C2conn[i] -= mean_x2 * mean_x2;
    }

    C1_data[t] = std::move(C1);
    C2raw_data[t] = std::move(C2raw);
    C2conn_data[t] = std::move(C2conn);
    C3_data[t] = std::move(C3);

    instanton_counts[t] = count_zero_crossings(path);
  }

  // ===============================
  // Export per-trial correlators (for jackknife in Python)
  // ===============================
  {
    const std::string trials_file = cooled ? "data/fig6_trials_correlators.csv"
                                           : "data/fig4_trials_correlators.csv";

    std::ofstream out_trials(trials_file);
    out_trials << "trial,tau,C1,C2raw,C2conn,C3\n";

    // Save each trial separately to allow external resampling analysis.
    for (int t = 0; t < trials; ++t) {
      for (int i = 0; i < N; ++i) {
        out_trials << t << "," << i * a << "," << C1_data[t][i] << ","
                   << C2raw_data[t][i] << "," << C2conn_data[t][i] << ","
                   << C3_data[t][i] << "\n";
      }
    }
  }

  // ===============================
  // Save Fig. 3 histogram data
  // ===============================
  if (!cooled) {
    std::ofstream out_pos(output_prefix + "_positions.csv");
    out_pos << "x\n";

    // Export all sampled positions for the probability histogram P(x).
    for (double xi : all_positions) {
      out_pos << xi << "\n";
    }
  }

  // ===============================
  // Ensemble averages: Fig.4 (uncooled) and Fig.6 (cooled)
  // ===============================
  {
    std::vector<double> C1_mean(N, 0.0), C1_err(N, 0.0);
    std::vector<double> C2raw_mean(N, 0.0), C2raw_err(N, 0.0);
    std::vector<double> C2conn_mean(N, 0.0), C2conn_err(N, 0.0);
    std::vector<double> C3_mean(N, 0.0), C3_err(N, 0.0);

    // Average each correlator independently at every Euclidean time.
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

      double var1 = 0.0, var2raw = 0.0, var2conn = 0.0, var3 = 0.0;

      // Estimate statistical errors from trial-to-trial fluctuations.
      for (int t = 0; t < trials; ++t) {
        var1 += std::pow(C1_data[t][i] - C1_mean[i], 2);
        var2raw += std::pow(C2raw_data[t][i] - C2raw_mean[i], 2);
        var2conn += std::pow(C2conn_data[t][i] - C2conn_mean[i], 2);
        var3 += std::pow(C3_data[t][i] - C3_mean[i], 2);
      }

      if (trials > 1) {
        var1 /= (trials - 1);
        var2raw /= (trials - 1);
        var2conn /= (trials - 1);
        var3 /= (trials - 1);
      }

      C1_err[i] = std::sqrt(var1 / trials);
      C2raw_err[i] = std::sqrt(var2raw / trials);
      C2conn_err[i] = std::sqrt(var2conn / trials);
      C3_err[i] = std::sqrt(var3 / trials);
    }

    const std::string filename = cooled ? "data/fig6_cooled_correlators.csv"
                                        : "data/fig4_quantum_correlators.csv";

    std::ofstream out(filename);
    out << "tau,C1,C1_err,C2raw,C2raw_err,C2conn,C2conn_err,C3,C3_err\n";

    // Export ensemble-averaged correlators with standard errors.
    for (int i = 0; i < N; ++i) {
      out << i * a << "," << C1_mean[i] << "," << C1_err[i] << ","
          << C2raw_mean[i] << "," << C2raw_err[i] << "," << C2conn_mean[i]
          << "," << C2conn_err[i] << "," << C3_mean[i] << "," << C3_err[i]
          << "\n";
    }
  }

  // ===============================
  // Instanton statistics
  // ===============================
  std::vector<double> inst_vals(instanton_counts.begin(),
                                instanton_counts.end());

  const double inst_mean = mean(inst_vals);
  const double inst_stderr = standard_error(inst_vals);

  const double inst_density = inst_mean / beta;

  std::cout << "[✓] Ensemble (" << trials << " configs, cooled=" << cooled
            << ")  avg instantons = " << inst_mean << " ± " << inst_stderr
            << "  → density = " << inst_density << "\n";
}