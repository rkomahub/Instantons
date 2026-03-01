#include "cooling_evolution.hpp"
#include "instanton.hpp"
#include "metropolis.hpp"
#include <fstream>
#include <vector>

/**
 * @file cooling_evolution.cpp
 * @brief Study of instanton density under cooling.
 *
 * Cooling removes short-distance quantum fluctuations by
 * accepting only updates that decrease the Euclidean action.
 *
 * As cooling proceeds:
 *   - Quantum noise disappears rapidly (short time scale)
 *   - Instanton–anti-instanton annihilation occurs more slowly
 *
 * Monitoring the number of zero crossings as a function of
 * cooling sweeps reveals the separation of these scales.
 */

void run_cooling_evolution(const Lattice &original, int max_sweeps, double a,
                           const std::string &output_filename) {
  std::ofstream out(output_filename);
  out << "n_cool,n_inst,density\n";

  // Work on a copy so the original configuration is unchanged.
  Lattice cooled = original;

  // Metropolis object used in "cooling mode"
  Metropolis evolver(cooled);

  for (int n_cool = 1; n_cool <= max_sweeps; ++n_cool) {

    // Perform exactly one cooling sweep.
    evolver.cool(1);

    // Count tunneling events via zero crossings.
    const int n_inst = count_zero_crossings(cooled.get_path());

    // Total Euclidean time extent.
    const double beta = cooled.size() * a;

    // Instanton density n_{I+A} = N_{I+A} / beta.
    const double density = static_cast<double>(n_inst) / beta;

    out << n_cool << "," << n_inst << "," << density << "\n";
  }
}