#include "metropolis.hpp"
#include "parameters.hpp"
#include "potential.hpp"
#include <cmath>
#include <random>

/**
 * @file metropolis.cpp
 * @brief Implementation of Metropolis updates for the Euclidean lattice action.
 *
 * The target distribution is
 *
 *     P[x] ∝ exp(−S_E[x])
 *
 * with discretized Euclidean action
 *
 *     S_E = Σ_i [ (x_i − x_{i−1})^2 / (4a) + a V(x_i) ]
 *
 * where V(x) = (x^2 − eta^2)^2 and periodic boundary conditions are assumed.
 */

Metropolis::Metropolis(Lattice &lattice, std::mt19937 &gen)
    : x(lattice), gen(gen) {}

/**
 * @brief Perform one Metropolis sweep over all lattice sites.
 *
 * For each site i, the algorithm proposes
 *
 *     x_i → x_i + δx ,   δx ~ N(0, dx_width)
 *
 * and accepts with probability
 *
 *     min(1, exp(−ΔS))
 *
 * where ΔS is the change in the *local* Euclidean action terms that depend on
 * x_i.
 */
void Metropolis::step() {
  std::normal_distribution<double> dx_dist(0.0, params::dx_width);
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

  for (int i = 0; i < x.size(); ++i) {
    const int ip = (i + 1) % x.size();
    const int im = (i - 1 + x.size()) % x.size();

    const double x_old = x[i];

    // Local action contribution involving x_i:
    //   (x_i − x_{i−1})^2/(4a) + (x_{i+1} − x_i)^2/(4a) + a V(x_i)
    const double S_old =
        (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
            (4.0 * params::a) +
        params::a * potential(x[i], params::eta);

    x[i] += dx_dist(gen);

    const double S_new =
        (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
            (4.0 * params::a) +
        params::a * potential(x[i], params::eta);

    const double dS = S_new - S_old;

    // If u > exp(−ΔS), then reject; else accept.
    if (prob_dist(gen) > std::exp(-dS)) {
      x[i] = x_old;
    }
  }
}

/**
 * @brief Perform n_sweeps of cooling updates.
 *
 * Cooling accepts a proposal only if it lowers the local Euclidean action:
 * if ΔS < 0 then accept, else reject.
 */
void Metropolis::cool(int n_sweeps) {
  std::normal_distribution<double> dx_dist(0.0, params::dx_width_cool);

  for (int sweep = 0; sweep < n_sweeps; ++sweep) {
    for (int i = 0; i < x.size(); ++i) {
      const int ip = (i + 1) % x.size();
      const int im = (i - 1 + x.size()) % x.size();

      const double x_old = x[i];

      const double S_old =
          (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
              (4.0 * params::a) +
          params::a * potential(x[i], params::eta);

      x[i] += dx_dist(gen);

      const double S_new =
          (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
              (4.0 * params::a) +
          params::a * potential(x[i], params::eta);

      const double dS = S_new - S_old;

      // Cooling acceptance rule
      if (dS > 0.0) {
        x[i] = x_old;
      }
    }
  }
}