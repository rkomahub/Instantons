#include "analysis_driver.hpp"
#include "cooling_evolution.hpp"
#include "ensemble.hpp"
#include "fig14_ia_interaction.hpp"
#include "heating.hpp"
#include "iilm.hpp"
#include "instanton.hpp"
#include "io.hpp"
#include "lattice.hpp"
#include "metropolis.hpp"
#include "observables.hpp"
#include "parameters.hpp"
#include "potential.hpp"
#include "qmidens.hpp"
#include "rilm.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

struct Configuration {
  std::string label;
  std::vector<double> path;

  Configuration(std::string l, std::vector<double> p)
      : label(std::move(l)), path(std::move(p)) {}
};

void run_basic_analysis(std::mt19937 &gen) {
  Lattice lattice(params::N, params::eta, params::hot_start, gen);
  Metropolis evolver(lattice, gen);

  for (int sweep = 0; sweep < params::sweeps; ++sweep) {
    evolver.step();
  }

  // --- Store configurations ---
  std::vector<Configuration> configs;

  // Quantum (uncooled) path
  configs.emplace_back("quantum", lattice.get_path());

  // Cooled path (copy, then apply cooling)
  Lattice cooled = lattice;
  Metropolis cooled_evolver(cooled, gen);
  cooled_evolver.cool(200); // 200 cooling sweeps
  configs.emplace_back("cooled", cooled.get_path());

  // --- Analyze and save each configuration ---
  for (const auto &config : configs) {
    std::string base = "data/" + config.label;

    // Save the path
    save_path_to_csv(config.path, base + "_path.csv", params::a);
    // Save the correlator
    auto corr = compute_correlator(config.path);
    save_correlator_to_csv(corr, params::a, base + "_correlator.csv");

    // Instanton count
    int n_inst = count_zero_crossings(config.path);
    double beta = config.path.size() * params::a;
    double density = static_cast<double>(n_inst) / beta;

    std::cout << "[✓] " << config.label << ": " << n_inst
              << " instantons, density = " << density << "\n";
  }
}

void run_cooling_evolution_analysis(std::mt19937 &gen) {
  Lattice lattice(params::N, params::eta, params::hot_start, gen);
  Metropolis evolver(lattice, gen);
  for (int sweep = 0; sweep < params::sweeps; ++sweep) {
    evolver.step();
  }

  run_cooling_evolution(lattice, 200, params::a,
                        "data/instanton_density_vs_ncool.csv", gen);
}

void run_ensemble_analysis(std::mt19937 &gen) {
  int trials = 50;
  run_ensemble_average(trials, false, "data/ensemble_quantum", gen);
  run_ensemble_average(trials, true, "data/ensemble_cooled", gen);
}

void run_fig7_analysis(std::mt19937 &gen) {
  const std::vector<double> etas = {1.4, 1.5, 1.6};

  const int Nconf = 50;
  const int skip = 200;
  const int ncool_max = 200;

  std::ofstream out("data/fig7_long.csv");
  out << "eta,conf,n_cool,n_inst,density,action,s_per_inst\n";

  for (double eta_val : etas) {
    params::eta = eta_val;

    Lattice lat(params::N, params::eta, params::hot_start, gen);
    Metropolis evo(lat, gen);

    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    for (int k = 0; k < Nconf; ++k) {
      for (int s = 0; s < skip; ++s)
        evo.step();

      Lattice cooled = lat;
      Metropolis cool_evo(cooled, gen);

      const double beta = cooled.size() * params::a;

      for (int n_cool = 0; n_cool <= ncool_max; ++n_cool) {
        if (n_cool > 0)
          cool_evo.cool(1);

        const int n_inst = count_zero_crossings(cooled.get_path());
        const double density = double(n_inst) / beta;

        const double S =
            compute_action(cooled.get_path(), params::a, params::eta);
        const double s_per_inst = (n_inst > 0) ? (S / n_inst) : NAN;

        out << eta_val << "," << k << "," << n_cool << "," << n_inst << ","
            << density << "," << S << "," << s_per_inst << "\n";
      }
    }
  }
}

void run_fig8_analysis(const std::vector<double> &etas, std::mt19937 &gen) {
  // --- BLUE (MC + 10 cooling sweeps)
  const int Nconf = 100;
  const int skip = 300;
  const int ncool = 10;

  std::ofstream outB("data/fig8_cooling10.csv");
  outB << "eta,density_mean,density_err\n";

  // --- RED (non-Gaussian switching)
  const int sweeps = 80;
  const double dx = params::dx_width;
  const int n_alpha_fine = 21;
  const int n_alpha_coarse = 11;

  std::ofstream outR("data/fig8_qmidens.csv");
  outR << "eta,density_ng_mean,density_ng_err\n";

  for (double eta_val : etas) {
    params::eta = eta_val;

    // ===== BLUE: full MC -> cool 10 -> count
    Lattice lat(params::N, params::eta, params::hot_start, gen);
    Metropolis evo(lat, gen);
    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    const double beta = params::N * params::a;

    std::vector<double> dens;
    dens.reserve(Nconf);

    for (int k = 0; k < Nconf; ++k) {
      for (int s = 0; s < skip; ++s)
        evo.step();

      Lattice cooled = lat;
      Metropolis cool_evo(cooled, gen);
      cool_evo.cool(ncool);

      const int n_inst = count_zero_crossings(cooled.get_path());
      dens.push_back(double(n_inst) / beta);
    }

    double mean_blue = 0.0;
    for (double x : dens)
      mean_blue += x;
    mean_blue /= double(dens.size());

    double var_blue = 0.0;
    for (double x : dens)
      var_blue += (x - mean_blue) * (x - mean_blue);
    if (dens.size() > 1)
      var_blue /= double(dens.size() - 1);

    const double err_blue = std::sqrt(var_blue / double(dens.size()));
    outB << eta_val << "," << mean_blue << "," << err_blue << "\n";

    // ===== RED: non-Gaussian corrected density
    const int Nrep = 8; // increase if errors are noisy

    std::vector<double> vals;
    vals.reserve(Nrep);

    for (int r = 0; r < Nrep; ++r) {
      // independent RNG stream per replicate
      std::mt19937 gen_r(gen()); // draws a seed from the main generator

      const double n_ng = compute_qmidens_corrected_density(
          sweeps, dx, n_alpha_fine, n_alpha_coarse, gen_r);

      vals.push_back(n_ng);
    }

    double mean_ng = 0.0;
    for (double v : vals)
      mean_ng += v;
    mean_ng /= double(vals.size());

    double var_ng = 0.0;
    for (double v : vals)
      var_ng += (v - mean_ng) * (v - mean_ng);
    if (vals.size() > 1)
      var_ng /= double(vals.size() - 1);

    const double err_ng = std::sqrt(var_ng / double(vals.size()));

    outR << eta_val << "," << mean_ng << "," << err_ng << "\n";
  }
}

void run_fig9_analysis(std::mt19937 &gen) {
  std::cout << "[📊] Running Fig.9 switching path generation...\n";
  run_fig9_paths(gen);
}

void run_rilm_analysis(std::mt19937 &gen) {
  const int trials = 50;
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;

  int n_inst = 4; // keep fixed for now
  if (n_inst % 2 != 0)
    ++n_inst; // safety

  std::vector<std::vector<double>> C1_data(trials), C2raw_data(trials),
      C2conn_data(trials), C3_data(trials);

  for (int t = 0; t < trials; ++t) {
    auto path = generate_rilm_path(N, a, eta, n_inst, gen);

    if (t == 0) {
      auto [mn, mx] = std::minmax_element(path.begin(), path.end());
      std::cout << "trial0 path min=" << *mn << " max=" << *mx << "\n";
      save_path_to_csv(path, "data/fig12_rilm_path.csv", a);
    }

    const double mean_x2 = compute_moment(path, 2);

    auto C1 = compute_correlator_power(path, 1);
    auto C2raw = compute_correlator_power(path, 2);
    auto C3 = compute_correlator_power(path, 3);

    auto C2conn = C2raw;
    for (int i = 0; i < N; ++i)
      C2conn[i] -= mean_x2 * mean_x2;

    // store (now all vectors have size N)
    C1_data[t] = std::move(C1);
    C2raw_data[t] = std::move(C2raw);
    C2conn_data[t] = std::move(C2conn);
    C3_data[t] = std::move(C3);
  }

  // ------------------------------
  // Export per-trial file (AFTER loop)
  // ------------------------------
  std::ofstream outT("data/fig10_trials_correlators.csv");
  outT << "trial,tau,C1,C2raw,C2conn,C3\n";
  for (int t = 0; t < trials; ++t) {
    for (int i = 0; i < N; ++i) {
      outT << t << "," << i * a << "," << C1_data[t][i] << ","
           << C2raw_data[t][i] << "," << C2conn_data[t][i] << ","
           << C3_data[t][i] << "\n";
    }
  }

  // ------------------------------
  // Compute ensemble averages (AFTER loop)
  // ------------------------------
  std::vector<double> C1_mean(N, 0), C1_err(N, 0);
  std::vector<double> C2raw_mean(N, 0), C2raw_err(N, 0);
  std::vector<double> C2conn_mean(N, 0), C2conn_err(N, 0);
  std::vector<double> C3_mean(N, 0), C3_err(N, 0);

  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < trials; ++t) {
      C1_mean[i] += C1_data[t][i];
      C2raw_mean[i] += C2raw_data[t][i];
      C2conn_mean[i] += C2conn_data[t][i];
      C3_mean[i] += C3_data[t][i];
    }
    C1_mean[i] /= trials;
    C2raw_mean[i] /= trials;
    C2conn_mean[i] /= trials;
    C3_mean[i] /= trials;

    double v1 = 0, v2r = 0, v2c = 0, v3 = 0;
    for (int t = 0; t < trials; ++t) {
      v1 += std::pow(C1_data[t][i] - C1_mean[i], 2);
      v2r += std::pow(C2raw_data[t][i] - C2raw_mean[i], 2);
      v2c += std::pow(C2conn_data[t][i] - C2conn_mean[i], 2);
      v3 += std::pow(C3_data[t][i] - C3_mean[i], 2);
    }
    if (trials > 1) {
      v1 /= (trials - 1);
      v2r /= (trials - 1);
      v2c /= (trials - 1);
      v3 /= (trials - 1);
    }
    C1_err[i] = std::sqrt(v1 / trials);
    C2raw_err[i] = std::sqrt(v2r / trials);
    C2conn_err[i] = std::sqrt(v2c / trials);
    C3_err[i] = std::sqrt(v3 / trials);
  }

  std::ofstream out("data/fig10_rilm_correlators.csv");
  out << "tau,C1,C1_err,C2raw,C2raw_err,C2conn,C2conn_err,C3,C3_err\n";
  for (int i = 0; i < N; ++i) {
    out << i * a << "," << C1_mean[i] << "," << C1_err[i] << ","
        << C2raw_mean[i] << "," << C2raw_err[i] << "," << C2conn_mean[i] << ","
        << C2conn_err[i] << "," << C3_mean[i] << "," << C3_err[i] << "\n";
  }

  std::cout << "[✓] Fig.10 RILM ensemble generated\n";
}

void run_heated_rilm_analysis(std::mt19937 &gen) {
  const int trials = 100; // before 50
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;

  // --- RILM parameters (match your Fig.10 choice)
  int n_inst = 4;
  if (n_inst % 2 != 0)
    ++n_inst;

  // --- heating parameters (Fig.12/13)
  const int n_white_sweeps = 10;
  const double fluct_sigma = 0.2;

  const int n_heat_sweeps = 300;    // metro sweeps (before 10, 250)
  const double dx_width_heat = 0.3; // tune acceptance (before 0.05, 0.2)

  auto Vdd = [eta](double xc) { return 12.0 * xc * xc - 4.0 * eta * eta; };
  std::normal_distribution<double> dx_dist(0.0, dx_width_heat);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  // --- storage for Fig.13 correlators (heated ensemble)
  std::vector<std::vector<double>> C1_data(trials), C2raw_data(trials),
      C2conn_data(trials), C3_data(trials);

  long long n_prop_tot = 0, n_acc_tot = 0;

  for (int t = 0; t < trials; ++t) {

    // 1) fresh backbone per trial (this is what Fig.13 means)
    auto x_cl = generate_rilm_path(N, a, eta, n_inst, gen);

    // 2) for t==0 export Fig.12 typical paths (white + metro) from this
    // backbone
    if (t == 0) {
      save_path_to_csv(x_cl, "data/fig12_rilm_path.csv", a);

      // (A) white-noise heated typical path
      auto x_white = x_cl;
      for (int s = 0; s < n_white_sweeps; ++s) {
        apply_gaussian_fluctuations(x_white, fluct_sigma, gen);
      }
      save_path_to_csv(x_white, "data/fig12_heated_rilm_path.csv", a);

      // (B) metro heated typical path
      auto x_metro = x_cl;
      for (int sweep = 0; sweep < n_heat_sweeps; ++sweep) {
        for (int i = 0; i < N; ++i) {
          const int ip = (i + 1) % N;
          const int im = (i - 1 + N) % N;

          const double x_old = x_metro[i];
          const double kappa = Vdd(x_cl[i]);
          const double dx0 = x_old - x_cl[i];

          const double S_old = (std::pow(x_old - x_metro[im], 2) +
                                std::pow(x_metro[ip] - x_old, 2)) /
                                   (4.0 * a) +
                               0.5 * a * kappa * dx0 * dx0;

          const double x_new = x_old + dx_dist(gen);
          const double dx1 = x_new - x_cl[i];

          const double S_new = (std::pow(x_new - x_metro[im], 2) +
                                std::pow(x_metro[ip] - x_new, 2)) /
                                   (4.0 * a) +
                               0.5 * a * kappa * dx1 * dx1;

          const double dS = S_new - S_old;

          ++n_prop_tot;
          if (u01(gen) <= std::exp(-dS)) {
            x_metro[i] = x_new;
            ++n_acc_tot;
          }
        }
      }
      save_path_to_csv(x_metro, "data/fig12_heated_metro_rilm_path.csv", a);
    }

    // 3) metro-heated configuration for correlators (Fig.13 ensemble)
    std::vector<double> path = x_cl;
    for (int sweep = 0; sweep < n_heat_sweeps; ++sweep) {
      for (int i = 0; i < N; ++i) {
        const int ip = (i + 1) % N;
        const int im = (i - 1 + N) % N;

        const double x_old = path[i];
        const double kappa = Vdd(x_cl[i]);
        const double dx0 = x_old - x_cl[i];

        const double S_old =
            (std::pow(x_old - path[im], 2) + std::pow(path[ip] - x_old, 2)) /
                (4.0 * a) +
            0.5 * a * kappa * dx0 * dx0;

        const double x_new = x_old + dx_dist(gen);
        const double dx1 = x_new - x_cl[i];

        const double S_new =
            (std::pow(x_new - path[im], 2) + std::pow(path[ip] - x_new, 2)) /
                (4.0 * a) +
            0.5 * a * kappa * dx1 * dx1;

        const double dS = S_new - S_old;

        ++n_prop_tot;
        if (u01(gen) <= std::exp(-dS)) {
          path[i] = x_new;
          ++n_acc_tot;
        }
      }
    }

    // 4) correlators (same as Fig.10/4)
    const double mean_x2 = compute_moment(path, 2);

    auto C1 = compute_correlator_power(path, 1);
    auto C2raw = compute_correlator_power(path, 2);
    auto C3 = compute_correlator_power(path, 3);

    auto C2conn = C2raw;
    for (int i = 0; i < N; ++i)
      C2conn[i] -= mean_x2 * mean_x2;

    C1_data[t] = std::move(C1);
    C2raw_data[t] = std::move(C2raw);
    C2conn_data[t] = std::move(C2conn);
    C3_data[t] = std::move(C3);
  }

  // ---- export per-trial for jackknife (Fig.13)
  {
    std::ofstream outT("data/fig13_trials_correlators.csv");
    outT << "trial,tau,C1,C2raw,C2conn,C3\n";
    for (int t = 0; t < trials; ++t) {
      for (int i = 0; i < N; ++i) {
        outT << t << "," << i * a << "," << C1_data[t][i] << ","
             << C2raw_data[t][i] << "," << C2conn_data[t][i] << ","
             << C3_data[t][i] << "\n";
      }
    }
  }

  // ---- ensemble averages + stderr (Fig.13)
  std::vector<double> C1_mean(N, 0), C1_err(N, 0);
  std::vector<double> C2raw_mean(N, 0), C2raw_err(N, 0);
  std::vector<double> C2conn_mean(N, 0), C2conn_err(N, 0);
  std::vector<double> C3_mean(N, 0), C3_err(N, 0);

  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < trials; ++t) {
      C1_mean[i] += C1_data[t][i];
      C2raw_mean[i] += C2raw_data[t][i];
      C2conn_mean[i] += C2conn_data[t][i];
      C3_mean[i] += C3_data[t][i];
    }
    C1_mean[i] /= trials;
    C2raw_mean[i] /= trials;
    C2conn_mean[i] /= trials;
    C3_mean[i] /= trials;

    double v1 = 0, v2r = 0, v2c = 0, v3 = 0;
    for (int t = 0; t < trials; ++t) {
      v1 += std::pow(C1_data[t][i] - C1_mean[i], 2);
      v2r += std::pow(C2raw_data[t][i] - C2raw_mean[i], 2);
      v2c += std::pow(C2conn_data[t][i] - C2conn_mean[i], 2);
      v3 += std::pow(C3_data[t][i] - C3_mean[i], 2);
    }
    if (trials > 1) {
      v1 /= (trials - 1);
      v2r /= (trials - 1);
      v2c /= (trials - 1);
      v3 /= (trials - 1);
    }

    C1_err[i] = std::sqrt(v1 / trials);
    C2raw_err[i] = std::sqrt(v2r / trials);
    C2conn_err[i] = std::sqrt(v2c / trials);
    C3_err[i] = std::sqrt(v3 / trials);
  }

  {
    std::ofstream out("data/fig13_heated_rilm_correlators.csv");
    out << "tau,C1,C1_err,C2raw,C2raw_err,C2conn,C2conn_err,C3,C3_err\n";
    for (int i = 0; i < N; ++i) {
      out << i * a << "," << C1_mean[i] << "," << C1_err[i] << ","
          << C2raw_mean[i] << "," << C2raw_err[i] << "," << C2conn_mean[i]
          << "," << C2conn_err[i] << "," << C3_mean[i] << "," << C3_err[i]
          << "\n";
    }
  }

  const double acc_rate =
      (n_prop_tot > 0) ? double(n_acc_tot) / double(n_prop_tot) : 0.0;

  std::cout << "[✓] Fig.12 paths written (t=0) and Fig.13 heated-RILM "
               "correlators written (trials="
            << trials << ", metro acceptance=" << acc_rate << ")\n";
}

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

void run_fig15_analysis(std::mt19937 & /*gen*/) {
  std::cout << "[📊] Running Fig.15 relaxed-IA family scan...\n";

  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;
  const double S0 = 4.0 * std::pow(eta, 3) / 3.0;

  const double tauI = 0.25 * beta;

  const std::vector<double> targets = {2.0, 1.8, 1.6, 1.4, 1.2, 1.0,
                                       0.8, 0.6, 0.4, 0.2, 0.1};

  std::ofstream out("data/fig15_paths.csv");
  out << "ratio,tau,x,s\n";

  // Scan many initial IA separations.
  const double tauIA_min = 0.10;
  const double tauIA_max = 4.00;
  const double dtau_scan = 0.02;

  // For each initial separation, apply a moderate amount of flow.
  const int nflow = 400;
  const double eps = 0.002;

  // Store best representative for each target.
  std::vector<std::vector<double>> best_paths(targets.size());
  std::vector<double> best_ratio(targets.size(), 0.0);
  std::vector<double> best_diff(targets.size(), 1e100);
  std::vector<bool> found(targets.size(), false);

  for (double tauIA0 = tauIA_min; tauIA0 <= tauIA_max + 1e-12;
       tauIA0 += dtau_scan) {
    const double tauA0 = std::fmod(tauI + tauIA0, beta);

    auto x = build_sum_ansatz_ia(N, a, eta, tauI, tauA0);

    for (int step = 0; step <= nflow; ++step) {
      const double S = compute_action(x, a, eta);
      const double ratio = S / S0;

      for (size_t k = 0; k < targets.size(); ++k) {
        const double diff = std::fabs(ratio - targets[k]);
        if (diff < best_diff[k]) {
          best_diff[k] = diff;
          best_ratio[k] = ratio;
          best_paths[k] = x;
          found[k] = true;
        }
      }

      gradient_flow_step(x, a, eta, eps);
    }
  }

  // Write one representative curve per target.
  for (size_t k = 0; k < targets.size(); ++k) {
    if (!found[k] || best_paths[k].empty()) {
      std::cout << "[warn] no representative found for target " << targets[k]
                << "\n";
      continue;
    }

    const auto &x = best_paths[k];

    std::cout << "target " << targets[k]
              << "  best actual ratio = " << best_ratio[k]
              << "  |diff| = " << best_diff[k] << "\n";

    for (int i = 0; i < N; ++i) {
      const int im = (i - 1 + N) % N;
      const double dx = (x[i] - x[im]) / a;
      const double s = 0.25 * dx * dx + potential(x[i], eta);
      const double tau = i * a;

      out << targets[k] << "," << tau << "," << x[i] << "," << s << "\n";
    }
  }

  std::cout << "[✓] data/fig15_paths.csv written\n";
}

void run_fig16_analysis(std::mt19937 &gen) { run_fig14_analysis(gen); }

void run_iilm_analysis(std::mt19937 &gen) {
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;

  const int n_inst = 10;
  const double tau_core = 0.3;

  // --- generate collective coordinates (reproducible)
  auto cfg = generate_iilm_config(n_inst, beta, tau_core, gen);

  // --- build lattice path from cfg
  auto iilm_path = build_iilm_path(N, a, eta, cfg);

  save_path_to_csv(iilm_path, "data/iilm_path.csv", a);
  save_correlator_to_csv(compute_correlator(iilm_path), a,
                         "data/iilm_correlator.csv");

  const int counted = count_zero_crossings(iilm_path);
  const double density = static_cast<double>(counted) / beta;

  std::cout << "[✓] IILM: placed " << n_inst << ", counted " << counted
            << ", density = " << density << "\n";
}

void run_fig17_analysis(std::mt19937 &gen) {
  std::cout << "[📊] Running Fig.17 IILM trajectory of instanton positions...\n";

  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;
  const double beta = N * a;

  const int n_inst_total =
      20; // should be >= 20 to have first 10 + and first 10 -
  const double tau_core = 0.3; // same as your run_iilm_analysis
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

void run_cooling_eta_scan(const std::vector<double> &etas, std::mt19937 &gen) {
  std::ofstream out("data/cooling_eta_scan.csv");
  out << "eta,density\n";
  for (double eta : etas) {
    Lattice lat(params::N, eta, /*hot_start=*/true, gen);
    Metropolis evo(lat, gen);
    for (int s = 0; s < params::sweeps; ++s)
      evo.step();

    // do ~10 cooling sweeps
    Metropolis cool(lat, gen);
    for (int c = 0; c < 10; ++c)
      cool.cool(1);

    int n = count_zero_crossings(lat.get_path());
    double beta = params::N * params::a;
    out << eta << "," << (double(n) / beta) << "\n";
  }
}
