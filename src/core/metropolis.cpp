#include "core/metropolis.hpp"
#include "core/potential.hpp"
#include "utils/parameters.hpp"
#include "utils/periodic.hpp"

#include <cmath>
#include <random>

namespace {

// Compute the part of the action affected by a local update at site i.
double local_action(const Lattice &x, int i) {
  const int N = x.size();

  const int ip = periodic_next(i, N);
  const int im = periodic_prev(i, N);

  const double dx_left = x[i] - x[im];
  const double dx_right = x[ip] - x[i];

  return (dx_left * dx_left + dx_right * dx_right) / (4.0 * params::a) +
         params::a * potential(x[i], params::eta);
}

} // namespace

// Bind the Metropolis updater to an existing lattice and RNG.
Metropolis::Metropolis(Lattice &lattice, std::mt19937 &gen)
    : x(lattice), gen(gen) {}

// Perform one full Metropolis sweep over all lattice sites.
void Metropolis::step() {
  std::normal_distribution<double> dx_dist(0.0, params::dx_width);
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  const int N = x.size();

  for (int i = 0; i < N; ++i) {

    const double x_old = x[i];

    // --- OLD action ---
    const double S_old = local_action(x, i);

    // --- proposal ---
    x[i] += dx_dist(gen);

    // --- NEW action ---
    const double S_new = local_action(x, i);

    const double dS = S_new - S_old;

    // --- Metropolis accept/reject ---
    if (prob_dist(gen) > std::exp(-dS)) {
      x[i] = x_old;
    }
  }
}

// Perform cooling sweeps, accepting only action-lowering updates.
void Metropolis::cool(int n_sweeps) {
  std::normal_distribution<double> dx_dist(0.0, params::dx_width_cool);

  for (int sweep = 0; sweep < n_sweeps; ++sweep) {
    for (int i = 0; i < x.size(); ++i) {
      const int ip = (i + 1) % x.size();
      const int im = (i - 1 + x.size()) % x.size();

      const double x_old = x[i];

      const double S_old = local_action(x, i);

      // Propose a small local deformation of the cooled path.
      x[i] += dx_dist(gen);

      const double S_new = local_action(x, i);

      const double dS = S_new - S_old;

      // Reject any update that increases the Euclidean action.
      if (dS > 0.0) {
        x[i] = x_old;
      }
    }
  }
}