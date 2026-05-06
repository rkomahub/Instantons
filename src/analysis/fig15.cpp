#include "analysis/fig15.hpp"

#include "core/observables.hpp"
#include "core/potential.hpp"
#include "models/fig14_ia_interaction.hpp"
#include "utils/parameters.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

void run_fig15_analysis(std::mt19937 & /*gen*/) {
  std::cout << "[📊] Running Fig.15 relaxed-IA family scan...\n";

  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;
  const double S0 = 4.0 * std::pow(eta, 3) / 3.0;

  const double tauI = 0.25 * beta;

  const std::vector<double> targets = {2.0, 1.8, 1.6, 1.4, 1.2, 1.0,
                                       0.8, 0.6, 0.4, 0.2, 0.1};

  std::ofstream out("data/fig15_paths.csv");
  out << "ratio,tau,x,s\n";

  // Scan many initial IA separations.
  const double tauIA_min = 0.10;
  const double tauIA_max = 4.00;
  const double dtau_scan = 0.02;

  // For each initial separation, apply a moderate amount of flow.
  const int nflow = 400;
  const double eps = 0.002;

  // Store best representative for each target.
  std::vector<std::vector<double>> best_paths(targets.size());
  std::vector<double> best_ratio(targets.size(), 0.0);
  std::vector<double> best_diff(targets.size(), 1e100);
  std::vector<bool> found(targets.size(), false);

  for (double tauIA0 = tauIA_min; tauIA0 <= tauIA_max + 1e-12;
       tauIA0 += dtau_scan) {
    const double tauA0 = std::fmod(tauI + tauIA0, beta);

    auto x = build_sum_ansatz_ia(N, a, eta, tauI, tauA0);

    for (int step = 0; step <= nflow; ++step) {
      const double S = compute_action(x, a, eta);
      const double ratio = S / S0;

      for (size_t k = 0; k < targets.size(); ++k) {
        const double diff = std::fabs(ratio - targets[k]);
        if (diff < best_diff[k]) {
          best_diff[k] = diff;
          best_ratio[k] = ratio;
          best_paths[k] = x;
          found[k] = true;
        }
      }

      gradient_flow_step(x, a, eta, eps);
    }
  }

  // Write one representative curve per target.
  for (size_t k = 0; k < targets.size(); ++k) {
    if (!found[k] || best_paths[k].empty()) {
      std::cout << "[warn] no representative found for target " << targets[k]
                << "\n";
      continue;
    }

    const auto &x = best_paths[k];

    std::cout << "target " << targets[k]
              << "  best actual ratio = " << best_ratio[k]
              << "  |diff| = " << best_diff[k] << "\n";

    for (int i = 0; i < N; ++i) {
      const int im = (i - 1 + N) % N;
      const double dx = (x[i] - x[im]) / a;
      const double s = 0.25 * dx * dx + potential(x[i], eta);
      const double tau = i * a;

      out << targets[k] << "," << tau << "," << x[i] << "," << s << "\n";
    }
  }

  std::cout << "[✓] data/fig15_paths.csv written\n";
}