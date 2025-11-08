#include "qmidens.hpp"
#include "instanton.hpp"
#include "io.hpp"
#include "observables.hpp"
#include "parameters.hpp"
#include "potential.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

namespace {

double S_gaussian(const std::vector<double> &x, double a, double eta) {
  double omega = 4.0 * eta;
  double Sg = 0.0;
  int N = x.size();
  for (int i = 0; i < N; ++i) {
    int ip = (i + 1) % N;
    double dx = x[ip] - x[i];
    double Vg = 0.5 * omega * omega * x[i] * x[i];
    Sg += dx * dx / (4.0 * a) + a * Vg;
  }
  return Sg;
}

double S_full(const std::vector<double> &x, double a, double eta) {
  double S = 0.0;
  int N = x.size();
  for (int i = 0; i < N; ++i) {
    int ip = (i + 1) % N;
    double dx = x[ip] - x[i];
    double Vf = std::pow(x[i] * x[i] - eta * eta, 2);
    S += dx * dx / (4.0 * a) + a * Vf;
  }
  return S;
}

double S_alpha(const std::vector<double> &x, double a, double eta,
               double alpha) {
  double Sg = S_gaussian(x, a, eta);
  double Sf = S_full(x, a, eta);
  return (1 - alpha) * Sg + alpha * Sf;
}

void metropolis_alpha(std::vector<double> &x, double a, double eta,
                      double alpha, int sweeps, double dx_width) {
  std::mt19937 gen(std::random_device{}());
  std::normal_distribution<double> dx_dist(0.0, dx_width);
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  int N = x.size();

  for (int sweep = 0; sweep < sweeps; ++sweep) {
    for (int i = 0; i < N; ++i) {
      int ip = (i + 1) % N;
      int im = (i - 1 + N) % N;

      double x_old = x[i];
      double dS_old = S_alpha(x, a, eta, alpha);

      x[i] += dx_dist(gen);
      double dS_new = S_alpha(x, a, eta, alpha);

      double dS = dS_new - dS_old;
      if (prob_dist(gen) > std::exp(-dS)) {
        x[i] = x_old;
      }
    }
  }
}

// Generate a random path constrained to n zero-crossings
std::vector<double> generate_constrained_path(int N, double eta,
                                              int n_zero_crossings) {
  std::vector<double> x(N);
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> dist(-eta, eta);

  do {
    for (int i = 0; i < N; ++i) {
      x[i] = dist(gen);
    }
  } while (count_zero_crossings(x) != n_zero_crossings);

  return x;
}

// Compute <ΔS> at fixed α in sector with `n_inst` zero crossings
double delta_S_alpha(int n_inst, double alpha, int sweeps, double dx_width) {
  int N = params::N;
  double a = params::a;
  double eta = params::eta;

  auto path = generate_constrained_path(N, eta, n_inst);
  metropolis_alpha(path, a, eta, alpha, sweeps, dx_width);

  double deltaS = S_full(path, a, eta) - S_gaussian(path, a, eta);
  return deltaS;
}

} // end anonymous namespace

double simpson_integral(int n_alpha, int sweeps, double dx,
                        std::ofstream *log = nullptr) {
  if (n_alpha % 2 == 0)
    ++n_alpha; // Simpson’s rule requires odd

  double h = 1.0 / (n_alpha - 1);
  double result = 0.0;

  for (int i = 0; i < n_alpha; ++i) {
    double alpha = static_cast<double>(i) / (n_alpha - 1);
    double ds1 = delta_S_alpha(1, alpha, sweeps, dx);
    double ds0 = delta_S_alpha(0, alpha, sweeps, dx);
    double diff = ds1 - ds0;

    if (log) {
      *log << alpha << "," << ds1 << "," << ds0 << "," << diff << "\n";
    }

    if (i == 0 || i == n_alpha - 1)
      result += diff;
    else if (i % 2 == 1)
      result += 4.0 * diff;
    else
      result += 2.0 * diff;
  }

  return result * h / 3.0;
}

void run_qmidens_analysis() {
  std::cout << "[📊] Running adiabatic switching (non-Gaussian correction)...\n";

  int sweeps = 1000;
  double dx = params::dx_width;

  int n_alpha_fine = 21;   // h = 1 / 20
  int n_alpha_coarse = 11; // 2h = 1 / 10

  // Optional: log the fine grid values to CSV
  std::ofstream out("data/qmidens_integrand.csv");
  out << "alpha,DeltaS1,DeltaS0,diff\n";

  double I_fine = simpson_integral(n_alpha_fine, sweeps, dx, &out);
  double I_coarse = simpson_integral(n_alpha_coarse, sweeps, dx);

  double deltaS_richardson =
      (16.0 * I_fine - I_coarse) / 15.0; // Richardson extrapolation

  const double S0 = 4.0 * std::pow(params::eta, 3) / 3.0;
  const double pref = 8.0 * std::pow(params::eta, 2.5) * std::sqrt(2.0 / M_PI);
  double n_gauss = pref * std::exp(-S0);
  double n_corrected = n_gauss * std::exp(-deltaS_richardson);

  std::cout << "[✓] I_fine   = " << I_fine << "\n";
  std::cout << "[✓] I_coarse = " << I_coarse << "\n";
  std::cout << "[✓] Richardson extrapolated ΔS = " << deltaS_richardson << "\n";
  std::cout << "[✓] n_gauss = " << n_gauss << ", n_corrected = " << n_corrected
            << "\n";
}

double compute_qmidens_corrected_density(int sweeps, double dx_width,
                                         int n_alpha_fine, int n_alpha_coarse) {
  // Ensure odd Simpson counts
  if (n_alpha_fine % 2 == 0)
    ++n_alpha_fine;
  if (n_alpha_coarse % 2 == 0)
    ++n_alpha_coarse;

  // Fine + coarse Simpson, then Richardson
  double I_fine = simpson_integral(n_alpha_fine, sweeps, dx_width);
  double I_coarse = simpson_integral(n_alpha_coarse, sweeps, dx_width);
  double deltaS_richardson = (16.0 * I_fine - I_coarse) / 15.0;

  // Semiclassical (Gaussian) prefactor and classical action:
  // S0 = 4 η^3 / 3  and  pref = 8 η^(5/2) * sqrt(2/π)
  const double S0 = 4.0 * std::pow(params::eta, 3) / 3.0;
  const double pref = 8.0 * std::pow(params::eta, 2.5) * std::sqrt(2.0 / M_PI);

  const double n_gauss = pref * std::exp(-S0);
  const double n_corrected = n_gauss * std::exp(-deltaS_richardson);
  return n_corrected;
}