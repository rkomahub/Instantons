#include "analysis/fig8.hpp"

#include "analysis/qmidens.hpp"
#include "core/instanton.hpp"
#include "core/lattice.hpp"
#include "core/metropolis.hpp"
#include "utils/parameters.hpp"
#include "utils/statistics.hpp"

#include <fstream>
#include <random>
#include <vector>

void run_fig8_analysis(const std::vector<double> &etas, std::mt19937 &gen) {
  // --- BLUE (MC + 10 cooling sweeps)
  const int Nconf = 100;
  const int skip = 300;
  const int ncool = 10;

  std::ofstream outB("data/fig8_cooling10.csv");
  outB << "eta,density_mean,density_err\n";

  // --- RED (non-Gaussian switching)
  const int sweeps = 80;
  const double dx = params::dx_width;
  const int n_alpha_fine = 21;
  const int n_alpha_coarse = 11;

  std::ofstream outR("data/fig8_qmidens.csv");
  outR << "eta,density_ng_mean,density_ng_err\n";

  for (double eta_val : etas) {
    params::eta = eta_val;

    // ===== BLUE: full MC -> cool 10 -> count
    Lattice lat(params::N, params::eta, params::hot_start, gen);
    Metropolis evo(lat, gen);
    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    const double beta = params::N * params::a;

    std::vector<double> dens;
    dens.reserve(Nconf);

    for (int k = 0; k < Nconf; ++k) {
      for (int s = 0; s < skip; ++s)
        evo.step();

      Lattice cooled = lat;
      Metropolis cool_evo(cooled, gen);
      cool_evo.cool(ncool);

      const int n_inst = count_zero_crossings(cooled.get_path());
      dens.push_back(double(n_inst) / beta);
    }

    const double mean_blue = mean(dens);
    const double err_blue = standard_error(dens);
    outB << eta_val << "," << mean_blue << "," << err_blue << "\n";

    // ===== RED: non-Gaussian corrected density
    const int Nrep = 8; // increase if errors are noisy

    std::vector<double> vals;
    vals.reserve(Nrep);

    for (int r = 0; r < Nrep; ++r) {
      // independent RNG stream per replicate
      std::mt19937 gen_r(gen()); // draws a seed from the main generator

      const double n_ng = compute_qmidens_corrected_density(
          sweeps, dx, n_alpha_fine, n_alpha_coarse, gen_r);

      vals.push_back(n_ng);
    }

    const double mean_ng = mean(vals);
    const double err_ng = standard_error(vals);
    outR << eta_val << "," << mean_ng << "," << err_ng << "\n";
  }
}