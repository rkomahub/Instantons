#include "analysis/analysis_driver.hpp"
#include "analysis/cooling_evolution.hpp"
#include "analysis/ensemble.hpp"
#include "analysis/qmidens.hpp"
#include "core/instanton.hpp"
#include "core/lattice.hpp"
#include "core/metropolis.hpp"
#include "core/observables.hpp"
#include "core/potential.hpp"
#include "models/fig14_ia_interaction.hpp"
#include "models/heating.hpp"
#include "models/iilm.hpp"
#include "models/rilm.hpp"
#include "utils/io.hpp"
#include "utils/parameters.hpp"
#include "utils/statistics.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Small container used to label and store one Euclidean path.
struct Configuration {
  std::string label;
  std::vector<double> path;

  Configuration(std::string l, std::vector<double> p)
      : label(std::move(l)), path(std::move(p)) {}
};

// Run one full Metropolis simulation and compare quantum and cooled paths.
void run_basic_analysis(std::mt19937 &gen) {
  Lattice lattice(params::N, params::eta, params::hot_start, gen);
  Metropolis evolver(lattice, gen);

  // Equilibrate and sample the Euclidean path with Metropolis sweeps.
  for (int sweep = 0; sweep < params::sweeps; ++sweep) {
    evolver.step();
  }

  // --- Store configurations ---
  std::vector<Configuration> configs;

  // Quantum (uncooled) path
  configs.emplace_back("quantum", lattice.get_path());

  // Cooled path (copy, then apply cooling)
  Lattice cooled = lattice;
  Metropolis cooled_evolver(cooled, gen);
  cooled_evolver.cool(200); // 200 cooling sweeps
  configs.emplace_back("cooled", cooled.get_path());

  // --- Analyze and save each configuration ---
  for (const auto &config : configs) {
    std::string base = "data/" + config.label;

    // Save the path
    save_path_to_csv(config.path, base + "_path.csv", params::a);
    // Save the correlator
    auto corr = compute_correlator(config.path);
    save_correlator_to_csv(corr, params::a, base + "_correlator.csv");

    // Estimate instanton content through zero crossings.
    int n_inst = count_zero_crossings(config.path);
    double beta = config.path.size() * params::a;
    double density = static_cast<double>(n_inst) / beta;

    std::cout << "[✓] " << config.label << ": " << n_inst
              << " instantons, density = " << density << "\n";
  }
}

// Measure how the instanton density changes under repeated cooling sweeps.
void run_cooling_evolution_analysis(std::mt19937 &gen) {
  Lattice lattice(params::N, params::eta, params::hot_start, gen);
  Metropolis evolver(lattice, gen);

  // Generate the initial quantum configuration before cooling.
  for (int sweep = 0; sweep < params::sweeps; ++sweep) {
    evolver.step();
  }

  run_cooling_evolution(lattice, 200, params::a,
                        "data/instanton_density_vs_ncool.csv", gen);
}

// Compute ensemble averages for quantum and cooled configurations.
void run_ensemble_analysis(std::mt19937 &gen) {
  int trials = 50;
  run_ensemble_average(trials, false, "data/ensemble_quantum", gen);
  run_ensemble_average(trials, true, "data/ensemble_cooled", gen);
}

// Generate and analyze an interacting instanton liquid configuration.
void run_iilm_analysis(std::mt19937 &gen) {
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;

  const int n_inst = 10;
  const double tau_core = 0.3;

  // --- generate collective coordinates (reproducible)
  auto cfg = generate_iilm_config(n_inst, beta, tau_core, gen);

  // --- build lattice path from cfg
  auto iilm_path = build_iilm_path(N, a, eta, cfg);

  // Save the generated IILM path and its Euclidean correlator.
  save_path_to_csv(iilm_path, "data/iilm_path.csv", a);
  save_correlator_to_csv(compute_correlator(iilm_path), a,
                         "data/iilm_correlator.csv");

  // Compare inserted instantons with the zero-crossing estimate.
  const int counted = count_zero_crossings(iilm_path);
  const double density = static_cast<double>(counted) / beta;

  std::cout << "[✓] IILM: placed " << n_inst << ", counted " << counted
            << ", density = " << density << "\n";
}

// Scan eta and measure the cooled instanton density after fixed cooling.
void run_cooling_eta_scan(const std::vector<double> &etas, std::mt19937 &gen) {
  std::ofstream out("data/cooling_eta_scan.csv");
  out << "eta,density\n";
  for (double eta : etas) {
    Lattice lat(params::N, eta, /*hot_start=*/true, gen);
    Metropolis evo(lat, gen);

    // Generate a quantum path for the current double-well parameter eta.
    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    // do ~10 cooling sweeps
    Metropolis cool(lat, gen);
    for (int c = 0; c < 10; ++c)
      cool.cool(1);

    // Store the instanton density extracted after cooling.
    int n = count_zero_crossings(lat.get_path());
    double beta = params::N * params::a;
    out << eta << "," << (double(n) / beta) << "\n";
  }
}