#include "iilm.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// -----------------------------
// Small periodic helpers
// -----------------------------
static inline double wrap_0_beta(double t, double beta) {
  t = std::fmod(t, beta);
  if (t < 0.0)
    t += beta;
  return t;
}

static inline double periodic_distance(double t1, double t2, double beta) {
  const double dt = std::fabs(t1 - t2);
  return std::min(dt, beta - dt);
}

static inline double nearest_image_dt(double dt, double beta) {
  // map dt to (-beta/2, beta/2]
  while (dt <= -0.5 * beta)
    dt += beta;
  while (dt > +0.5 * beta)
    dt -= beta;
  return dt;
}

// -----------------------------
// Core constraint
// -----------------------------
bool iilm_core_ok(const IILMConfig &cfg, double beta, double tau_core) {
  const int n = static_cast<int>(cfg.tau.size());
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (periodic_distance(cfg.tau[i], cfg.tau[j], beta) < tau_core) {
        return false;
      }
    }
  }
  return true;
}

// -----------------------------
// Initial configuration
// -----------------------------
IILMConfig generate_iilm_config(int n_inst, double beta, double tau_core,
                                std::mt19937 &gen) {
  if (n_inst <= 0)
    throw std::runtime_error("generate_iilm_config: n_inst <= 0");
  if (beta <= 0.0)
    throw std::runtime_error("generate_iilm_config: beta <= 0");
  if (tau_core < 0.0)
    throw std::runtime_error("generate_iilm_config: tau_core < 0");

  std::uniform_real_distribution<double> dist_tau(0.0, beta);
  std::uniform_int_distribution<int> dist_sign(0, 1);

  IILMConfig cfg;
  cfg.tau.reserve(n_inst);
  cfg.Q.reserve(n_inst);

  int attempts = 0;
  const int max_attempts = 200 * n_inst;

  while (static_cast<int>(cfg.tau.size()) < n_inst && attempts < max_attempts) {
    const double t = dist_tau(gen);

    bool too_close = false;
    for (double existing : cfg.tau) {
      if (periodic_distance(t, existing, beta) < tau_core) {
        too_close = true;
        break;
      }
    }

    if (!too_close) {
      cfg.tau.push_back(t);
      cfg.Q.push_back(dist_sign(gen) == 0 ? -1 : +1);
    }

    ++attempts;
  }

  if (static_cast<int>(cfg.tau.size()) != n_inst) {
    throw std::runtime_error("generate_iilm_config: failed to place all "
                             "objects (tau_core too large?)");
  }

  // Keep them in increasing tau for nicer bookkeeping/plotting (optional)
  // We sort by tau while carrying Q with it.
  std::vector<int> idx(n_inst);
  for (int i = 0; i < n_inst; ++i)
    idx[i] = i;
  std::sort(idx.begin(), idx.end(),
            [&](int a, int b) { return cfg.tau[a] < cfg.tau[b]; });

  IILMConfig sorted;
  sorted.tau.resize(n_inst);
  sorted.Q.resize(n_inst);
  for (int k = 0; k < n_inst; ++k) {
    sorted.tau[k] = cfg.tau[idx[k]];
    sorted.Q[k] = cfg.Q[idx[k]];
  }
  return sorted;
}

// -----------------------------
// Build lattice path from cfg
// -----------------------------
std::vector<double> build_iilm_path(int N, double a, double eta,
                                    const IILMConfig &cfg) {
  if (N <= 0)
    throw std::runtime_error("build_iilm_path: N <= 0");
  if (a <= 0.0)
    throw std::runtime_error("build_iilm_path: a <= 0");
  if (cfg.tau.size() != cfg.Q.size())
    throw std::runtime_error("build_iilm_path: tau/Q size mismatch");

  const double beta = N * a;
  const double omega = 4.0 * eta;

  std::vector<double> x(N, 0.0);

  for (int i = 0; i < N; ++i) {
    const double tau = i * a;

    double sum = 0.0;
    for (size_t j = 0; j < cfg.tau.size(); ++j) {
      const double dt = nearest_image_dt(tau - cfg.tau[j], beta);
      sum += double(cfg.Q[j]) * std::tanh(0.5 * omega * dt);
    }

    x[i] = eta * (sum - 1.0);
  }

  return x;
}

// -----------------------------
// Proposal move for Markov chain
// -----------------------------
bool propose_move_tau(IILMConfig &cfg, double beta, double tau_core,
                      double step, std::mt19937 &gen) {
  if (cfg.tau.empty())
    return false;
  if (step <= 0.0)
    throw std::runtime_error("propose_move_tau: step <= 0");

  std::uniform_int_distribution<int> pick(0,
                                          static_cast<int>(cfg.tau.size()) - 1);
  std::uniform_real_distribution<double> dist_d(-step, step);

  const int j = pick(gen);
  const double old = cfg.tau[j];
  const double trial = wrap_0_beta(old + dist_d(gen), beta);

  // core check only against others
  for (size_t k = 0; k < cfg.tau.size(); ++k) {
    if (static_cast<int>(k) == j)
      continue;
    if (periodic_distance(trial, cfg.tau[k], beta) < tau_core) {
      return false; // immediately reject proposal
    }
  }

  cfg.tau[j] = trial;
  return true;
}