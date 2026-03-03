#include "ensemble.hpp"
#include "instanton.hpp"
#include "lattice.hpp"
#include "metropolis.hpp"
#include "observables.hpp"
#include "parameters.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

/**
 * @file ensemble.cpp
 * @brief Ensemble averaging of correlators and instanton observables.
 *
 * The function runs multiple independent Monte Carlo simulations,
 * optionally applies cooling, and estimates expectation values by:
 *
 *   mean( C(τ) ) = (1/T) Σ_t C_t(τ)
 *
 * together with the standard error of the mean:
 *
 *   stderr( C(τ) ) = sqrt( Var(C(τ)) / T )
 *
 * where T = trials.
 */

void run_ensemble_average(int trials, bool cooled,
                          const std::string &output_prefix, std::mt19937 &gen) {

  const int N = params::N;
  const double a = params::a;
  const double beta = N * a;

  // For Fig. 3
  std::vector<double> all_positions;
  if (!cooled) {
    all_positions.reserve(static_cast<size_t>(trials) * static_cast<size_t>(N));
  }

  // For Fig. 4
  std::vector<std::vector<double>> C1_data(trials);
  std::vector<std::vector<double>> C2_data(trials);
  std::vector<std::vector<double>> C3_data(trials);

  std::vector<int> instanton_counts(trials);

  // ===============================
  // Generate ensemble
  // ===============================
  for (int t = 0; t < trials; ++t) {

    Lattice lattice(N, params::eta, /*hot_start=*/true, gen);
    Metropolis evolver(lattice, gen);

    for (int sweep = 0; sweep < params::sweeps; ++sweep) {
      evolver.step();
    }

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

    // --- Fig. 4 correlators ---
    auto C1 = compute_correlator_power(path, 1);
    auto C2 = compute_correlator_power(path, 2);
    auto C3 = compute_correlator_power(path, 3);

    // Connected subtraction for x^2
    double mean_x2 = compute_moment(path, 2);
    for (int i = 0; i < N; ++i) {
      C2[i] -= mean_x2 * mean_x2;
    }

    C1_data[t] = C1;
    C2_data[t] = C2;
    C3_data[t] = C3;

    instanton_counts[t] = count_zero_crossings(path);
  }

  // ===============================
  // Save Fig. 3 histogram data
  // ===============================
  if (!cooled) {
    std::ofstream out_pos(output_prefix + "_positions.csv");
    out_pos << "x\n";
    for (double xi : all_positions)
      out_pos << xi << "\n";
  }

  // ===============================
  // Fig. 4: ensemble averages
  // ===============================
  if (!cooled) {

    std::vector<double> C1_mean(N, 0.0), C1_err(N, 0.0);
    std::vector<double> C2_mean(N, 0.0), C2_err(N, 0.0);
    std::vector<double> C3_mean(N, 0.0), C3_err(N, 0.0);

    for (int i = 0; i < N; ++i) {

      for (int t = 0; t < trials; ++t) {
        C1_mean[i] += C1_data[t][i];
        C2_mean[i] += C2_data[t][i];
        C3_mean[i] += C3_data[t][i];
      }

      C1_mean[i] /= trials;
      C2_mean[i] /= trials;
      C3_mean[i] /= trials;

      double var1 = 0.0, var2 = 0.0, var3 = 0.0;

      for (int t = 0; t < trials; ++t) {
        var1 += std::pow(C1_data[t][i] - C1_mean[i], 2);
        var2 += std::pow(C2_data[t][i] - C2_mean[i], 2);
        var3 += std::pow(C3_data[t][i] - C3_mean[i], 2);
      }

      if (trials > 1) {
        var1 /= (trials - 1);
        var2 /= (trials - 1);
        var3 /= (trials - 1);
      }

      C1_err[i] = std::sqrt(var1 / trials);
      C2_err[i] = std::sqrt(var2 / trials);
      C3_err[i] = std::sqrt(var3 / trials);
    }

    std::ofstream out("data/fig4_quantum_correlators.csv");
    out << "tau,C1,C1_err,C2conn,C2conn_err,C3,C3_err\n";

    for (int i = 0; i < N; ++i) {
      out << i * a << "," << C1_mean[i] << "," << C1_err[i] << "," << C2_mean[i]
          << "," << C2_err[i] << "," << C3_mean[i] << "," << C3_err[i] << "\n";
    }
  }

  // ===============================
  // Instanton statistics
  // ===============================
  const double inst_mean =
      std::accumulate(instanton_counts.begin(), instanton_counts.end(), 0.0) /
      static_cast<double>(trials);

  const double inst_density = inst_mean / beta;

  double inst_var = 0.0;
  for (int t = 0; t < trials; ++t) {
    double d = instanton_counts[t] - inst_mean;
    inst_var += d * d;
  }

  if (trials > 1)
    inst_var /= (trials - 1);

  const double inst_stderr = std::sqrt(inst_var / trials);

  std::cout << "[✓] Ensemble (" << trials << " configs, cooled=" << cooled
            << ")  avg instantons = " << inst_mean << " ± " << inst_stderr
            << "  → density = " << inst_density << "\n";
}