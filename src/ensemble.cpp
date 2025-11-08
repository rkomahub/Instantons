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

void run_ensemble_average(int trials, bool cooled,
                          const std::string &output_prefix) {
  int N = params::N;
  double a = params::a;
  double beta = N * a;

  std::vector<std::vector<double>> correlator_data(trials);
  std::vector<int> instanton_counts(trials);

  for (int t = 0; t < trials; ++t) {
    Lattice lattice(N, params::eta, true); // hot start
    Metropolis evolver(lattice);

    for (int sweep = 0; sweep < params::sweeps; ++sweep)
      evolver.step();

    if (cooled) {
      Lattice cooled = lattice;
      Metropolis cooled_evolver(cooled);
      cooled_evolver.cool(200);
      lattice = cooled;
    }

    auto path = lattice.get_path();
    correlator_data[t] = compute_correlator(path);
    instanton_counts[t] = count_zero_crossings(path);
  }

  // Compute mean and stderr
  std::vector<double> mean(N, 0.0), stderr(N, 0.0);

  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < trials; ++t)
      mean[i] += correlator_data[t][i];
    mean[i] /= trials;

    for (int t = 0; t < trials; ++t)
      stderr[i] += std::pow(correlator_data[t][i] - mean[i], 2);
    stderr[i] = std::sqrt(stderr[i] / trials) / std::sqrt(trials);
  }

  // Save to CSV
  std::ofstream out_corr(output_prefix + "_correlator_avg.csv");
  out_corr << "tau,mean,std_err\n";
  for (int i = 0; i < N; ++i)
    out_corr << i * a << "," << mean[i] << "," << stderr[i] << "\n";

  // Instanton stats
  double inst_mean =
      std::accumulate(instanton_counts.begin(), instanton_counts.end(), 0.0) /
      trials;
  double inst_density = inst_mean / beta;

  double sumsq = 0.0;
  for (int t = 0; t < trials; ++t)
    sumsq += std::pow(correlator_data[t][i] - mean[i], 2);

  double var = sumsq / std::max(1, trials - 1); // unbiased sample variance
  stderr[i] = std::sqrt(var / trials);          // standard error of the mean

  std::cout << "[✓] Ensemble (" << trials << " configs, cooled=" << cooled
            << ")  avg instantons = " << inst_mean
            << " → density = " << inst_density << "\n";
}
