#include "analysis/fig14_16.hpp"

#include "core/instanton.hpp"
#include "core/lattice.hpp"
#include "core/metropolis.hpp"
#include "core/observables.hpp"
#include "models/fig14_ia_interaction.hpp"
#include "models/rilm.hpp"
#include "utils/parameters.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

void run_fig14_analysis(std::mt19937 &gen) {
  std::cout << "[📊] Running Fig.14 + Fig.16 combined analysis...\n";

  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;
  const double S0 = 4.0 * std::pow(eta, 3) / 3.0;

  // ============================================================
  // Part A: Fig.14 sum ansatz + streamline
  // ============================================================
  std::ofstream out_tauIA("data/fig14_sum_tauIA.csv");
  out_tauIA << "tauIA,Sint_over_S0\n";

  std::ofstream out_tauZ("data/fig14_sum_tauZ.csv");
  out_tauZ << "tauZ,Sint_over_S0\n";

  std::ofstream out_stream("data/fig14_streamline_tauZ.csv");
  out_stream << "tauZ,Sint_over_S0,step\n";

  const double tauI = 0.25 * beta;

  for (double tauIA = 0.05; tauIA <= 4.50; tauIA += 0.025) {
    const double tauA = std::fmod(tauI + tauIA, beta);

    auto x = build_sum_ansatz_ia(N, a, eta, tauI, tauA);

    const double S = compute_action(x, a, eta);
    const double Sint_over_S0 = S / S0 - 2.0;

    out_tauIA << tauIA << "," << Sint_over_S0 << "\n";

    const double tauZ = zero_crossing_separation(x, a);
    if (tauZ > 0.0) {
      out_tauZ << tauZ << "," << Sint_over_S0 << "\n";
    }
  }

  {
    const double eps = 0.002;
    const int nflow = 1200;

    for (double tauIA0 = 4.0; tauIA0 >= 0.10; tauIA0 -= 0.05) {
      const double tauA0 = std::fmod(tauI + tauIA0, beta);
      auto x = build_sum_ansatz_ia(N, a, eta, tauI, tauA0);

      for (int step = 0; step < nflow; ++step) {
        gradient_flow_step(x, a, eta, eps);
      }

      const double tauZ = zero_crossing_separation(x, a);
      if (tauZ < 0.0)
        continue;

      const double S = compute_action(x, a, eta);
      out_stream << tauZ << "," << (S / S0 - 2.0) << "," << nflow << "\n";
    }
  }

  // ============================================================
  // Part B: Fig.16 separation histogram + Fig.14 MC points
  // ============================================================
  const int trials_mc = 300;
  const int ncool = 10;

  const int trials_ref = 300;
  const int n_inst_ref = 4; // simple random reference ensemble

  const int nbins = 25;
  const double tau_min = 0.0;
  const double tau_max = beta / 2.0;
  const double binw = (tau_max - tau_min) / nbins;

  std::vector<int> hist_mc(nbins, 0);
  std::vector<int> hist_ref(nbins, 0);

  std::ofstream out_sep("data/fig16_separations.csv");
  out_sep << "source,trial,sep\n"; // source = mc or ref

  auto fill_hist = [&](double sep, std::vector<int> &hist) {
    if (sep < tau_min || sep >= tau_max)
      return;
    const int b = static_cast<int>((sep - tau_min) / binw);
    if (b >= 0 && b < nbins)
      hist[b]++;
  };

  // ---------------- MC cooled ensemble ----------------
  for (int t = 0; t < trials_mc; ++t) {
    Lattice lat(N, eta, params::hot_start, gen);
    Metropolis evo(lat, gen);

    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    Metropolis cool(lat, gen);
    for (int c = 0; c < ncool; ++c)
      cool.cool(1);

    const auto &path = lat.get_path();
    auto tau = zero_crossings_interp(path, a);

    if (tau.size() < 2)
      continue;

    for (size_t i = 0; i < tau.size(); ++i) {
      const double t1 = tau[i];
      const double t2 = tau[(i + 1) % tau.size()];

      const double sep =
          (i + 1 < tau.size()) ? (t2 - t1) : (beta - t1 + tau[0]);

      out_sep << "mc," << t << "," << sep << "\n";
      fill_hist(sep, hist_mc);
    }
  }

  // ---------------- random reference ensemble ----------------
  for (int t = 0; t < trials_ref; ++t) {
    auto path = generate_rilm_path(N, a, eta, n_inst_ref, gen);
    auto tau = zero_crossings_interp(path, a);

    if (tau.size() < 2)
      continue;

    for (size_t i = 0; i < tau.size(); ++i) {
      const double t1 = tau[i];
      const double t2 = tau[(i + 1) % tau.size()];

      const double sep =
          (i + 1 < tau.size()) ? (t2 - t1) : (beta - t1 + tau[0]);

      out_sep << "ref," << t << "," << sep << "\n";
      fill_hist(sep, hist_ref);
    }
  }

  // Save histogram file for Fig.16 plotting
  std::ofstream out_hist("data/fig16_histograms.csv");
  out_hist << "tau_center,n_mc,n_ref\n";

  // Save effective interaction points for Fig.14 MC scatter
  std::ofstream out_mc("data/fig14_mc_tauZ.csv");
  out_mc << "tauZ,Sint_over_S0\n";

  double total_mc = 0.0;
  double total_ref = 0.0;
  for (int b = 0; b < nbins; ++b) {
    total_mc += hist_mc[b];
    total_ref += hist_ref[b];
  }

  for (int b = 0; b < nbins; ++b) {
    const double tau_c = tau_min + (b + 0.5) * binw;
    const double n_mc = hist_mc[b];
    const double n_ref = hist_ref[b];

    out_hist << tau_c << "," << n_mc << "," << n_ref << "\n";

    if (n_mc > 0.0 && n_ref > 0.0 && total_mc > 0.0 && total_ref > 0.0) {
      const double p_mc = n_mc / total_mc;
      const double p_ref = n_ref / total_ref;

      const double Sint_eff = -std::log(p_mc / p_ref);
      const double Sint_over_S0 = Sint_eff / S0;

      out_mc << tau_c << "," << Sint_over_S0 << "\n";
    }
  }

  std::cout << "[✅] Wrote:\n"
            << "  data/fig14_sum_tauIA.csv\n"
            << "  data/fig14_sum_tauZ.csv\n"
            << "  data/fig14_streamline_tauZ.csv\n"
            << "  data/fig14_mc_tauZ.csv\n"
            << "  data/fig16_separations.csv\n"
            << "  data/fig16_histograms.csv\n";
}

void run_fig16_analysis(std::mt19937 &gen) { run_fig14_analysis(gen); }