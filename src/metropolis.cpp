#include "metropolis.hpp"
#include "parameters.hpp"
#include "potential.hpp"
#include <cmath>
#include <random>

Metropolis::Metropolis(Lattice &lattice)
    : x(lattice), gen(std::random_device{}()) {}

void Metropolis::step() {
  std::normal_distribution<double> dx_dist(0.0, params::dx_width);
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

  for (int i = 0; i < x.size(); ++i) {
    int ip = (i + 1) % x.size();
    int im = (i - 1 + x.size()) % x.size();

    double x_old = x[i];
    double S_old = (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
                       (4 * params::a) +
                   params::a * potential(x[i], params::eta);

    x[i] += dx_dist(gen);

    double S_new = (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
                       (4 * params::a) +
                   params::a * potential(x[i], params::eta);

    double dS = S_new - S_old;
    if (prob_dist(gen) > std::exp(-dS)) {
      x[i] = x_old;
    }
  }
}
