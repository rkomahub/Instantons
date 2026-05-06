#include "analysis/cooling_evolution.hpp"
#include "core/instanton.hpp"
#include "core/metropolis.hpp"
#include "core/observables.hpp"
#include "utils/parameters.hpp"

#include <cmath>
#include <fstream>
#include <random>

void run_cooling_evolution(const Lattice &original, int max_sweeps, double a,
                           const std::string &output_filename,
                           std::mt19937 &gen) {
  std::ofstream out(output_filename);
  out << "n_cool,n_inst,density,action,s_per_inst\n";

  // Work on a copy so the original configuration is unchanged.
  Lattice cooled = original;

  // Metropolis object used in "cooling mode"
  Metropolis evolver(cooled, gen);

  for (int n_cool = 1; n_cool <= max_sweeps; ++n_cool) {

    // Perform exactly one cooling sweep.
    evolver.cool(1);

    // Count tunneling events via zero crossings.
    const int n_inst = count_zero_crossings(cooled.get_path());

    // Total Euclidean time extent.
    const double beta = cooled.size() * a;

    // Instanton density n_{I+A} = N_{I+A} / beta.
    const double density = static_cast<double>(n_inst) / beta;

    double S = compute_action(cooled.get_path(), a, params::eta);

    double S_inst = NAN;

    if (n_inst > 0)
      S_inst = S / n_inst;

    out << n_cool << "," << n_inst << "," << density << "," << S << ","
        << S_inst << "\n";
  }
}