#include "lattice.hpp"
#include "metropolis.hpp"
#include "parameters.hpp"
#include "potential.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

double total_action(const std::vector<double> &x) {
  const int N = static_cast<int>(x.size());
  double S = 0.0;

  for (int i = 0; i < N; ++i) {
    int ip = (i + 1) % N;
    S += std::pow(x[ip] - x[i], 2) / (4.0 * params::a);
    S += params::a * potential(x[i], params::eta);
  }

  return S;
}

double local_action(const std::vector<double> &x, int i) {
  const int N = static_cast<int>(x.size());
  int ip = (i + 1) % N;
  int im = (i - 1 + N) % N;

  return (std::pow(x[i] - x[im], 2) + std::pow(x[ip] - x[i], 2)) /
             (4.0 * params::a) +
         params::a * potential(x[i], params::eta);
}

int main() {

  std::mt19937 gen(12345);

  // ---------- TEST 1: Locality ----------
  {
    const int N = 8;
    Lattice lat(N, params::eta, false, gen);

    for (int i = 0; i < N; ++i)
      lat[i] = 0.1 * i;

    auto x = lat.get_path();

    const int site = 3;
    const double dx = 0.01;

    // ΔS for a proposal at "site" with baseline configuration
    const double old_local = local_action(x, site);
    const double x_old = x[site];

    x[site] = x_old + dx;
    const double new_local = local_action(x, site);
    const double dS1 = new_local - old_local;

    // Restore baseline value at site
    x[site] = x_old;

    // Modify a distant site (should not affect ΔS at "site")
    x[6] += 10.0;

    const double old_local2 = local_action(x, site);

    x[site] = x_old + dx; // SAME proposal step
    const double new_local2 = local_action(x, site);
    const double dS2 = new_local2 - old_local2;

    assert(std::abs(dS1 - dS2) < 1e-12);
  }

  // ---------- TEST 2: Cooling monotonicity ----------
  {
    const int N = 32;
    Lattice lat(N, params::eta, true, gen);
    Metropolis evo(lat, gen);

    const double S_before = total_action(lat.get_path());

    evo.cool(1);

    const double S_after = total_action(lat.get_path());

    assert(S_after <= S_before + 1e-12);
  }

  std::cout << "[✓] All action tests passed.\n";
  return 0;
}