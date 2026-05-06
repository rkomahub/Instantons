#include "core/lattice.hpp"
#include "core/metropolis.hpp"
#include "core/potential.hpp"
#include "utils/parameters.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

static double local_action_site(const std::vector<double> &x, int i) {
  const int N = static_cast<int>(x.size());
  const int ip = (i + 1) % N;
  const int im = (i - 1 + N) % N;

  return (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
             (4.0 * params::a) +
         params::a * potential(x[i], params::eta);
}

int main() {

  // We use two generators with the same seed:
  // - gen_pred predicts the first proposal dx
  // - gen_run is passed to Metropolis::cool and must generate the same dx
  const unsigned seed = 424242;
  std::mt19937 gen_pred(seed);
  std::mt19937 gen_run(seed);

  const int N = 16;
  Lattice lat(N, params::eta, /*hot_start=*/false, gen_run);

  // Make the configuration deterministic (not relying on hot/cold start
  // randomness)
  for (int i = 0; i < N; ++i)
    lat[i] = 0.1 * i - 0.2;

  const int site = 0;

  // Predict the first cooling proposal at site=0
  std::normal_distribution<double> dx_dist(0.0, params::dx_width_cool);
  const double dx0 = dx_dist(gen_pred);

  auto x_before = lat.get_path();
  const double S_old = local_action_site(x_before, site);

  auto x_proposed = x_before;
  x_proposed[site] += dx0;
  const double S_new = local_action_site(x_proposed, site);

  const double dS = S_new - S_old;

  // Run exactly one cooling sweep (site loop starts from 0)
  Metropolis evo(lat, gen_run);
  evo.cool(1);

  const auto x_after = lat.get_path();

  if (dS > 0.0) {
    assert(std::abs(x_after[site] - x_before[site]) < 1e-15);
  } else {
    assert(std::abs(x_after[site] - x_proposed[site]) < 1e-15);
  }

  std::cout << "[✓] Cooling acceptance rule test passed.\n";
  return 0;
}