#include "analysis/analysis_driver.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {
// Classical single-instanton profile centered at tau0.
double x_instanton(double tau, double eta, double tau0) {
  const double u = 2.0 * eta * (tau - tau0);
  return eta * std::tanh(u);
}

// Gaussian fluctuation potential around the instanton path.
double V_gauss(double tau, double eta, double tau0) {
  const double u = 2.0 * eta * (tau - tau0);
  const double c = std::cosh(u);
  const double sech2 = 1.0 / (c * c);
  return 4.0 * eta * eta * (1.0 - 1.5 * sech2);
}
} // namespace

// Export the instanton profile and Gaussian effective potential for Fig. 11.
void run_fig11_analysis(double eta) {
  std::cout << "[📈] Running Fig. 11 export (Gaussian effective potential)...\n";

  std::filesystem::create_directories("data");
  const std::string path = "data/fig11_gaussian_effective_potential.csv";
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Cannot open output file: " + path);
  }

  const double tau0 = 0.0;

  const double T = 4.0; // window length (adjust freely for nicer plot)
  const int N = 800;    // resolution
  const double tau_min = -0.5 * T;
  const double tau_max = +0.5 * T;

  out << "tau,xI,Vg\n";

  // Sample both functions on a uniform Euclidean-time grid.
  for (int i = 0; i < N; ++i) {
    const double f = (N == 1) ? 0.0 : double(i) / double(N - 1);
    const double tau = tau_min + f * (tau_max - tau_min);
    out << tau << "," << x_instanton(tau, eta, tau0) << ","
        << V_gauss(tau, eta, tau0) << "\n";
  }

  std::cout << "  -> wrote " << path << "\n";
}