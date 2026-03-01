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
                          const std::string &output_prefix) {
  const int N = params::N;
  const double a = params::a;
  const double beta = N * a;

  std::vector<std::vector<double>> correlator_data(trials);
  std::vector<int> instanton_counts(trials);

  // Generate an ensemble of (approximately) independent configurations.
  for (int t = 0; t < trials; ++t) {
    Lattice lattice(N, params::eta, /*hot_start=*/true);
    Metropolis evolver(lattice);

    for (int sweep = 0; sweep < params::sweeps; ++sweep) {
      evolver.step();
    }

    // Optionally cool before measuring observables.
    if (cooled) {
      Lattice cooled_lat = lattice;
      Metropolis cooled_evolver(cooled_lat);
      cooled_evolver.cool(200);
      lattice = cooled_lat;
    }

    const auto &path = lattice.get_path();
    correlator_data[t] = compute_correlator(path);
    instanton_counts[t] = count_zero_crossings(path);
  }

  // Compute mean correlator and standard error of the mean for each τ.
  std::vector<double> mean(N, 0.0), stderr(N, 0.0);

  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < trials; ++t) {
      mean[i] += correlator_data[t][i];
    }
    mean[i] /= static_cast<double>(trials);

    double var = 0.0;
    for (int t = 0; t < trials; ++t) {
      const double d = correlator_data[t][i] - mean[i];
      var += d * d;
    }

    // If trials > 1, use unbiased sample variance; else var = 0.
    if (trials > 1) {
      var /= static_cast<double>(trials - 1);
    } else {
      var = 0.0;
    }

    stderr[i] = std::sqrt(var / static_cast<double>(trials));
  }

  // Save correlator statistics to CSV.
  std::ofstream out_corr(output_prefix + "_correlator_avg.csv");
  out_corr << "tau,mean,std_err\n";
  for (int i = 0; i < N; ++i) {
    out_corr << i * a << "," << mean[i] << "," << stderr[i] << "\n";
  }

  // Instanton density statistics from zero-crossing counts.
  const double inst_mean =
      std::accumulate(instanton_counts.begin(), instanton_counts.end(), 0.0) /
      static_cast<double>(trials);

  const double inst_density = inst_mean / beta;

  double inst_var = 0.0;
  for (int t = 0; t < trials; ++t) {
    const double d = static_cast<double>(instanton_counts[t]) - inst_mean;
    inst_var += d * d;
  }
  if (trials > 1) {
    inst_var /= static_cast<double>(trials - 1);
  } else {
    inst_var = 0.0;
  }

  const double inst_stderr = std::sqrt(inst_var / static_cast<double>(trials));

  std::cout << "[✓] Ensemble (" << trials << " configs, cooled=" << cooled
            << ")  avg instantons = " << inst_mean << " ± " << inst_stderr
            << "  → density = " << inst_density << "\n";
}