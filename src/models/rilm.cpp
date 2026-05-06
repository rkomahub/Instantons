#include "models/rilm.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

inline double min_image_dt(double dt, double beta) {
  dt -= beta * std::round(dt / beta);
  return dt;
}

// Return true if sorted centers satisfy nearest-neighbor hard-core on a circle
bool passes_hard_core(const std::vector<double> &tau_sorted, double beta,
                      double dmin) {
  const int n = static_cast<int>(tau_sorted.size());
  for (int j = 0; j < n - 1; ++j) {
    if ((tau_sorted[j + 1] - tau_sorted[j]) <= dmin) {
      return false;
    }
  }
  // circular neighbor: last -> first
  const double wrap_gap = beta - (tau_sorted[n - 1] - tau_sorted[0]);
  return (wrap_gap > dmin);
}

// Sample n centers uniformly and enforce hard-core via rejection
std::vector<double> sample_centers_hard_core(int n, double beta, double dmin,
                                             std::mt19937 &gen,
                                             int max_tries = 20000) {
  if (n <= 0) {
    throw std::runtime_error("sample_centers_hard_core: n must be > 0");
  }
  if (dmin <= 0.0) {
    throw std::runtime_error("sample_centers_hard_core: dmin must be > 0");
  }
  // simple necessary condition: total excluded length must fit in beta
  // (very rough, but useful to fail fast)
  if (n * dmin >= beta) {
    throw std::runtime_error(
        "sample_centers_hard_core: n*dmin >= beta (impossible constraints)");
  }

  std::uniform_real_distribution<double> dist_tau(0.0, beta);
  std::vector<double> tau(n);

  for (int attempt = 0; attempt < max_tries; ++attempt) {
    for (int j = 0; j < n; ++j)
      tau[j] = dist_tau(gen);
    std::sort(tau.begin(), tau.end());

    if (passes_hard_core(tau, beta, dmin)) {
      return tau;
    }
  }

  throw std::runtime_error(
      "sample_centers_hard_core: failed to sample centers (increase max_tries "
      "or reduce dmin)");
}

} // namespace

// =================================================
// product-ansatz hypothesis
// =================================================
std::vector<double> generate_rilm_path(int N, double a, double eta, int n_inst,
                                       std::mt19937 &gen) {
  if (n_inst <= 0)
    throw std::runtime_error("generate_rilm_path: n_inst must be > 0");
  if (n_inst % 2 != 0)
    throw std::runtime_error("generate_rilm_path: n_inst must be even");

  const double beta = N * a;
  const double omega = 4.0 * eta;
  const double dmin = 5.0 / omega;

  // centers (sorted, hard-core on circle)
  const std::vector<double> tau_i =
      sample_centers_hard_core(n_inst, beta, dmin, gen);

  // pick initial vacuum sign
  std::uniform_int_distribution<int> pm(0, 1);
  const int s0 = (pm(gen) == 0) ? -1 : +1;

  std::vector<double> x(N, 0.0);

  for (int i = 0; i < N; ++i) {
    const double tau = i * a;

    // smooth sign field: product of flips
    double s = double(s0);
    for (int j = 0; j < n_inst; ++j) {
      const double dt = min_image_dt(tau - tau_i[j], beta);
      s *= -std::tanh(0.5 * omega * dt);
    }

    x[i] = eta * s;
  }

  return x;
}

// =================================================
// sum-ansatz hypothesis
// =================================================
/*
std::vector<double> generate_rilm_path(int N, double a, double eta, int n_inst,
                                       std::mt19937 &gen) {
  if (n_inst <= 0) {
    throw std::runtime_error("generate_rilm_path: n_inst must be > 0");
  }
  if (n_inst % 2 != 0) {
    throw std::runtime_error("generate_rilm_path: n_inst must be even");
  }

  const double beta = N * a;
  const double omega = 4.0 * eta;

  // Hard-core (in Euclidean time)
  const double dmin = 5.0 / omega;

  // 1) sample centers with hard core, sorted
  const std::vector<double> tau_i =
      sample_centers_hard_core(n_inst, beta, dmin, gen);

  // diagnostic: minimal neighbor separation on circle
  double dt_min = 1e300;
  for (int j = 0; j < n_inst; ++j) {
    const double t1 = tau_i[j];
    const double t2 = (j == n_inst - 1) ? (tau_i[0] + beta) : tau_i[j + 1];
    dt_min = std::min(dt_min, t2 - t1);
  }

  // 2) charges: n/2 instantons (+1) and n/2 anti-instantons (-1), shuffled ONCE
  std::vector<int> Q;
  Q.reserve(n_inst);
  for (int k = 0; k < n_inst / 2; ++k)
    Q.push_back(+1);
  for (int k = 0; k < n_inst / 2; ++k)
    Q.push_back(-1);
  std::shuffle(Q.begin(), Q.end(), gen);

  // one-time signature print
  static bool once = false;
  if (!once) {
    once = true;
    std::cout << "[RILM] charges shuffled with net zero; dmin=" << dmin
              << " n_inst=" << n_inst << "\n";
  }
  std::cout << "[RILM] dt_min=" << dt_min << "  omega*dt_min=" << omega * dt_min
            << "\n";

  // 3) build path with SUM ansatz (periodic min-image separation)
  std::vector<double> x(N, 0.0);

  for (int i = 0; i < N; ++i) {
    const double tau = i * a;

    double sum = 0.0;
    for (int j = 0; j < n_inst; ++j) {
      const double dt = min_image_dt(tau - tau_i[j], beta);
      sum += Q[j] * eta * std::tanh(0.5 * omega * dt);
    }

    x[i] = sum - eta;
  }

  return x;
}
*/