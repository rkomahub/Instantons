#include "cooling_evolution.hpp"
#include "instanton.hpp"
#include "metropolis.hpp"
#include <fstream>
#include <vector>

void run_cooling_evolution(const Lattice &original, int max_sweeps, double a,
                           const std::string &output_filename) {
  std::ofstream out(output_filename);
  out << "n_cool,n_inst,density\n";

  Lattice cooled = original;
  Metropolis evolver(cooled);

  for (int n_cool = 1; n_cool <= max_sweeps; ++n_cool) {
    evolver.cool(1); // perform 1 sweep at a time

    int n_inst = count_zero_crossings(cooled.get_path());
    double beta = cooled.size() * a;
    double density = static_cast<double>(n_inst) / beta;

    out << n_cool << "," << n_inst << "," << density << "\n";
  }
}
