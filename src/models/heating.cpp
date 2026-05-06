#include "models/heating.hpp"
#include <random>

// Add Gaussian quantum fluctuations to an existing classical path.
void apply_gaussian_fluctuations(std::vector<double> &path, double sigma,
                                 std::mt19937 &gen) {
  std::normal_distribution<double> gauss(0.0, sigma);

  for (auto &x : path) {

    // Independently fluctuate each lattice site.
    x += gauss(gen);
  }
}