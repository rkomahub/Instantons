#include "potential.hpp"
#include <cmath>

double potential(double x, double eta) {
  return std::pow(x * x - eta * eta, 2);
}
