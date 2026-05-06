#include "analysis/fig7.hpp"

#include "core/instanton.hpp"
#include "core/lattice.hpp"
#include "core/metropolis.hpp"
#include "core/observables.hpp"
#include "utils/parameters.hpp"

#include <cmath>
#include <fstream>
#include <random>
#include <vector>

void run_fig7_analysis(std::mt19937 &gen) {
  const std::vector<double> etas = {1.4, 1.5, 1.6};

  const int Nconf = 50;
  const int skip = 200;
  const int ncool_max = 200;

  std::ofstream out("data/fig7_long.csv");
  out << "eta,conf,n_cool,n_inst,density,action,s_per_inst\n";

  for (double eta_val : etas) {
    params::eta = eta_val;

    Lattice lat(params::N, params::eta, params::hot_start, gen);
    Metropolis evo(lat, gen);

    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    for (int k = 0; k < Nconf; ++k) {
      for (int s = 0; s < skip; ++s)
        evo.step();

      Lattice cooled = lat;
      Metropolis cool_evo(cooled, gen);

      const double beta = cooled.size() * params::a;

      for (int n_cool = 0; n_cool <= ncool_max; ++n_cool) {
        if (n_cool > 0)
          cool_evo.cool(1);

        const int n_inst = count_zero_crossings(cooled.get_path());
        const double density = double(n_inst) / beta;

        const double S =
            compute_action(cooled.get_path(), params::a, params::eta);
        const double s_per_inst = (n_inst > 0) ? (S / n_inst) : NAN;

        out << eta_val << "," << k << "," << n_cool << "," << n_inst << ","
            << density << "," << S << "," << s_per_inst << "\n";
      }
    }
  }
}
