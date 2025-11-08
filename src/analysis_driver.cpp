#include "analysis_driver.hpp"
#include "cooling_evolution.hpp"
#include "ensemble.hpp"
#include "heating.hpp"
#include "iilm.hpp"
#include "instanton.hpp"
#include "io.hpp"
#include "lattice.hpp"
#include "metropolis.hpp"
#include "observables.hpp"
#include "parameters.hpp"
#include "qmidens.hpp"
#include "rilm.hpp"

#include <iostream>
#include <string>
#include <vector>

struct Configuration {
  std::string label;
  std::vector<double> path;

  Configuration(std::string l, std::vector<double> p)
      : label(std::move(l)), path(std::move(p)) {}
};

void run_basic_analysis() {
  Lattice lattice(params::N, params::eta, params::hot_start);
  Metropolis evolver(lattice);

  for (int sweep = 0; sweep < params::sweeps; ++sweep) {
    evolver.step();
  }

  // --- Store configurations ---
  std::vector<Configuration> configs;

  // Quantum (uncooked) path
  configs.emplace_back("quantum", lattice.get_path());

  // Cooled path (copy, then apply cooling)
  Lattice cooled = lattice;
  Metropolis cooled_evolver(cooled);
  cooled_evolver.cool(200);
  configs.emplace_back("cooled", cooled.get_path());

  // --- Analyze and save each configuration ---
  for (const auto &config : configs) {
    std::string base = "data/" + config.label;

    // Save the path
    save_path_to_csv(config.path, base + "_path.csv", params::a);
    // Save the correlator
    auto corr = compute_correlator(config.path);
    save_correlator_to_csv(corr, params::a, base + "_correlator.csv");

    // Instanton count
    int n_inst = count_zero_crossings(config.path);
    double beta = config.path.size() * params::a;
    double density = static_cast<double>(n_inst) / beta;

    std::cout << "[✓] " << config.label << ": " << n_inst
              << " instantons, density = " << density << "\n";
  }
}

void run_cooling_evolution_analysis() {
  Lattice lattice(params::N, params::eta, params::hot_start);
  Metropolis evolver(lattice);
  for (int sweep = 0; sweep < params::sweeps; ++sweep) {
    evolver.step();
  }

  run_cooling_evolution(lattice, 200, params::a,
                        "data/instanton_density_vs_ncool.csv");
}

void run_ensemble_analysis() {
  int trials = 50;
  run_ensemble_average(trials, false, "data/ensemble_quantum");
  run_ensemble_average(trials, true, "data/ensemble_cooled");
}

void run_rilm_analysis() {
  int N = params::N;
  double a = params::a;
  double eta = params::eta;
  int n_inst = 10;

  auto rilm_path = generate_rilm_path(N, a, eta, n_inst);

  save_path_to_csv(rilm_path, "data/rilm_path.csv", a);

  auto rilm_corr = compute_correlator(rilm_path);
  save_correlator_to_csv(rilm_corr, a, "data/rilm_correlator.csv");

  int counted = count_zero_crossings(rilm_path);
  double density = static_cast<double>(counted) / (N * a);

  std::cout << "[✓] RILM: placed " << n_inst << ", counted " << counted
            << ", density = " << density << "\n";
}

void run_heated_rilm_analysis() {
  int N = params::N;
  double a = params::a;
  double eta = params::eta;
  int n_inst = 10;
  double fluct_sigma = 0.2; // Tune this to match MC-like noise

  // 1. Generate clean RILM path
  auto rilm_path = generate_rilm_path(N, a, eta, n_inst);

  // 2. Copy and add Gaussian noise
  auto heated_path = rilm_path;
  apply_gaussian_fluctuations(heated_path, fluct_sigma);

  // 3. Save both
  save_path_to_csv(heated_path, "data/rilm_heated_path.csv", a);
  save_correlator_to_csv(compute_correlator(heated_path), a,
                         "data/rilm_heated_correlator.csv");

  int counted = count_zero_crossings(heated_path);
  double density = static_cast<double>(counted) / (N * a);

  std::cout << "[✓] Heated RILM: placed " << n_inst << ", counted " << counted
            << ", density = " << density << "\n";
}

void run_iilm_analysis() {
  int N = params::N;
  double a = params::a;
  double eta = params::eta;
  int n_inst = 10;
  double tau_core = 0.3;

  auto iilm_path = generate_iilm_path(N, a, eta, n_inst, tau_core);

  save_path_to_csv(iilm_path, "data/iilm_path.csv", a);
  save_correlator_to_csv(compute_correlator(iilm_path), a,
                         "data/iilm_correlator.csv");

  int counted = count_zero_crossings(iilm_path);
  double density = static_cast<double>(counted) / (N * a);

  std::cout << "[✓] IILM: placed " << n_inst << ", counted " << counted
            << ", density = " << density << "\n";
}

void run_qmidens_analysis() { ::run_qmidens_analysis(); }

void run_cooling_eta_scan(const std::vector<double> &etas) {
  std::ofstream out("/mnt/data/cooling_eta_scan.csv");
  out << "eta,density\n";
  for (double eta : etas) {
    Lattice lat(params::N, eta, /*hot_start=*/true);
    Metropolis evo(lat);
    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    // do ~10 cooling sweeps
    Metropolis cool(lat);
    for (int c = 0; c < 10; ++c)
      cool.cool(1);

    int n = count_zero_crossings(lat.get_path());
    double beta = params::N * params::a;
    out << eta << "," << (double(n) / beta) << "\n";
  }
}
