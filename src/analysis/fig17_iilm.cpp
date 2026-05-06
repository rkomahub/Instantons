#include "analysis/fig17_iilm.hpp"

#include "core/observables.hpp"
#include "models/iilm.hpp"
#include "utils/parameters.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

void run_fig17_analysis(std::mt19937 &gen) {
  std::cout << "[📊] Running Fig.17 IILM trajectory of instanton positions...\n";

  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;

  const int n_inst_total =
      20; // should be >= 20 to have first 10 + and first 10 -
  const double tau_core = 0.3; // same as run_iilm_analysis
  const int n_cfg = 3000;
  const double step_tau = 0.15; // proposal size in tau (tune acceptance)
  const int thin = 1;           // record every config; increase if needed

  // initial config + initial path/action
  IILMConfig cfg = generate_iilm_config(n_inst_total, beta, tau_core, gen);
  auto x = build_iilm_path(N, a, eta, cfg);
  double S = compute_action(x, a, eta);

  std::uniform_real_distribution<double> u01(0.0, 1.0);

  std::ofstream out("data/fig17_positions.csv");
  out << "conf,kind,k,tau\n"; // kind=+1 instanton, kind=-1 anti-instanton,
                              // k=0..9

  long long n_prop = 0, n_acc = 0;

  for (int conf = 0; conf < n_cfg; ++conf) {
    // --- Metropolis update in collective coordinates
    IILMConfig trial = cfg;

    ++n_prop;
    const bool proposed =
        propose_move_tau(trial, beta, tau_core, step_tau, gen);
    if (proposed) {
      auto x_trial = build_iilm_path(N, a, eta, trial);
      const double S_trial = compute_action(x_trial, a, eta);
      const double dS = S_trial - S;

      if (u01(gen) <= std::exp(-dS)) {
        cfg = std::move(trial);
        x = std::move(x_trial);
        S = S_trial;
        ++n_acc;
      }
    }

    // --- record first 10 I and first 10 A
    if ((conf % thin) == 0) {
      std::vector<double> tau_I, tau_A;
      tau_I.reserve(cfg.tau.size());
      tau_A.reserve(cfg.tau.size());

      for (size_t j = 0; j < cfg.tau.size(); ++j) {
        if (cfg.Q[j] == +1)
          tau_I.push_back(cfg.tau[j]);
        else
          tau_A.push_back(cfg.tau[j]);
      }

      std::sort(tau_I.begin(), tau_I.end());
      std::sort(tau_A.begin(), tau_A.end());

      const int nI = std::min<int>(10, (int)tau_I.size());
      const int nA = std::min<int>(10, (int)tau_A.size());

      for (int k = 0; k < nI; ++k)
        out << conf << "," << +1 << "," << k << "," << tau_I[k] << "\n";

      for (int k = 0; k < nA; ++k)
        out << conf << "," << -1 << "," << k << "," << tau_A[k] << "\n";
    }
  }

  const double acc = (n_prop > 0) ? double(n_acc) / double(n_prop) : 0.0;
  std::cout << "[✓] Wrote data/fig17_positions.csv, acceptance=" << acc << "\n";
}