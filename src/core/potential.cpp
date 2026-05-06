#include "core/potential.hpp"
#include <cmath>

// Double-well potential V(x) = (x^2 - eta^2)^2.
double potential(double x, double eta) {
  double y = x * x - eta * eta;
  return y * y;
}