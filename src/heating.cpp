#include "heating.hpp"
#include <random>

void apply_gaussian_fluctuations(std::vector<double> &path, double sigma) {
  std::mt19937 gen(std::random_device{}());
  std::normal_distribution<double> gauss(0.0, sigma);

  for (auto &x : path) {
    x += gauss(gen);
  }
}
