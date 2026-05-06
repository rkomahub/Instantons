#include "analysis/qmidens.hpp"
#include "core/instanton.hpp"
#include "core/observables.hpp"
#include "core/potential.hpp"
#include "utils/io.hpp"
#include "utils/parameters.hpp"
#include "utils/periodic.hpp"
#include "utils/statistics.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

namespace {

// ---------- potentials ----------
inline double V_alpha(double x, double eta, double alpha) {
  const double omega = 4.0 * eta;
  const double Vg = 0.25 * omega * omega * x * x;
  const double Vf = potential(x, eta);
  return (1.0 - alpha) * Vg + alpha * Vf;
}

inline double F_harmonic(double beta, double eta) {
  const double omega = 4.0 * eta;
  return (1.0 / beta) * std::log(2.0 * std::sinh(0.5 * beta * omega));
}

// ---------- actions ----------
double S_gaussian(const std::vector<double> &x, double a, double eta) {
  const double omega = 4.0 * eta;
  double Sg = 0.0;
  const int N = static_cast<int>(x.size());

  for (int i = 0; i < N; ++i) {
    const int ip = periodic_next(i, N);
    const double dx = x[ip] - x[i];
    const double Vg = 0.25 * omega * omega * x[i] * x[i];
    Sg += dx * dx / (4.0 * a) + a * Vg;
  }

  return Sg;
}

double S_full(const std::vector<double> &x, double a, double eta) {
  return compute_action(x, a, eta);
}

double S_alpha(const std::vector<double> &x, double a, double eta,
               double alpha) {
  const double Sg = S_gaussian(x, a, eta);
  const double Sf = S_full(x, a, eta);
  return (1.0 - alpha) * Sg + alpha * Sf;
}

// ---------- metropolis ----------
void metropolis_alpha(std::vector<double> &x, double a, double eta,
                      double alpha, int sweeps, double dx_width,
                      std::mt19937 &gen) {
  std::normal_distribution<double> dx_dist(0.0, dx_width);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  const int N = static_cast<int>(x.size());

  for (int sweep = 0; sweep < sweeps; ++sweep) {
    for (int i = 0; i < N; ++i) {
      const int ip = periodic_next(i, N);
      const int im = periodic_prev(i, N);

      const double x_old = x[i];
      const double x_new = x_old + dx_dist(gen);

      const double dx_left_new = x_new - x[im];
      const double dx_right_new = x[ip] - x_new;
      const double dx_left_old = x_old - x[im];
      const double dx_right_old = x[ip] - x_old;

      const double dS_kin =
          (dx_left_new * dx_left_new + dx_right_new * dx_right_new -
           dx_left_old * dx_left_old - dx_right_old * dx_right_old) /
          (4.0 * a);

      const double dS_pot =
          a * (V_alpha(x_new, eta, alpha) - V_alpha(x_old, eta, alpha));

      const double dS = dS_kin + dS_pot;

      if (u01(gen) < std::exp(-dS)) {
        x[i] = x_new;
      }
    }
  }
}

// ---------- constrained updates (used by Fig8) ----------
std::vector<double> init_vacuum_path(int N, double eta, int sign = -1) {
  return std::vector<double>(N, sign * eta);
}

std::vector<double> init_instanton_path(int N, double a, double eta, int i0) {
  std::vector<double> x(N);
  const double omega = 4.0 * eta;

  for (int i = 0; i < N; ++i) {
    const double tau = (i - i0) * a;
    x[i] = eta * std::tanh(0.5 * omega * tau);
  }

  x[i0] = 0.0;
  return x;
}

void metropolis_alpha_constrained(std::vector<double> &x, double a, double eta,
                                  double alpha, int sweeps, double dx_width,
                                  std::mt19937 &gen, int pinned_index,
                                  bool enforce_orientation) {
  std::normal_distribution<double> dx_dist(0.0, dx_width);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  const int N = static_cast<int>(x.size());
  const int i0 = pinned_index;

  for (int sweep = 0; sweep < sweeps; ++sweep) {
    for (int i = 0; i < N; ++i) {
      if (i0 >= 0 && i == i0) {
        x[i0] = 0.0;
        continue;
      }

      const int ip = periodic_next(i, N);
      const int im = periodic_prev(i, N);

      const double x_old = x[i];
      const double x_new = x_old + dx_dist(gen);

      const double dx_left_new = x_new - x[im];
      const double dx_right_new = x[ip] - x_new;
      const double dx_left_old = x_old - x[im];
      const double dx_right_old = x[ip] - x_old;

      const double dS_kin =
          (dx_left_new * dx_left_new + dx_right_new * dx_right_new -
           dx_left_old * dx_left_old - dx_right_old * dx_right_old) /
          (4.0 * a);

      const double dS_pot =
          a * (V_alpha(x_new, eta, alpha) - V_alpha(x_old, eta, alpha));

      const double dS = dS_kin + dS_pot;

      if (u01(gen) < std::exp(-dS)) {
        x[i] = x_new;

        if (i0 >= 0) {
          x[i0] = 0.0;
        }

        if (enforce_orientation && i0 >= 0) {
          const int im0 = periodic_prev(i0, N);
          const int ip0 = periodic_next(i0, N);

          if (!(x[im0] < 0.0 && x[ip0] > 0.0)) {
            x[i] = x_old;
            x[i0] = 0.0;
          }
        }
      }
    }
  }
}

// ---------- Fig5 helpers ----------

double measure_deltaS_on_path(std::vector<double> &path, double a, double eta,
                              double alpha, int n_therm, int n_meas, int stride,
                              double dx_width, std::mt19937 &gen) {
  metropolis_alpha(path, a, eta, alpha, n_therm, dx_width, gen);

  double sum = 0.0;
  int count = 0;

  for (int sweep = 0; sweep < n_meas; ++sweep) {
    metropolis_alpha(path, a, eta, alpha, 1, dx_width, gen);

    if (sweep % stride == 0) {
      const double deltaS = S_full(path, a, eta) - S_gaussian(path, a, eta);
      sum += deltaS;
      ++count;
    }
  }

  return (count > 0) ? sum / static_cast<double>(count) : 0.0;
}

double switching_integral_free_energy_beta(int n_alpha, int sweeps,
                                           double dx_width, int N_beta,
                                           std::mt19937 &gen,
                                           std::ofstream *log = nullptr) {
  if (n_alpha % 2 == 0) {
    ++n_alpha;
  }

  const double h = 1.0 / (n_alpha - 1);
  double result = 0.0;

  const double a = params::a;
  const double eta = params::eta;

  const int n_therm = sweeps;
  const int n_meas = 10 * sweeps;
  const int stride = 10;

  std::vector<double> path(N_beta, -eta);

  metropolis_alpha(path, a, eta, 0.0, 5 * sweeps, dx_width, gen);

  for (int i = 0; i < n_alpha; ++i) {
    const double alpha =
        static_cast<double>(i) / static_cast<double>(n_alpha - 1);

    const double ds = measure_deltaS_on_path(path, a, eta, alpha, n_therm,
                                             n_meas, stride, dx_width, gen);

    if (log) {
      *log << alpha << "," << ds << "\n";
    }

    if (i == 0 || i == n_alpha - 1) {
      result += ds;
    } else if (i % 2 == 1) {
      result += 4.0 * ds;
    } else {
      result += 2.0 * ds;
    }
  }

  return result * h / 3.0;
}

// ---------- sector generators ----------
std::vector<double> generate_constrained_path(int N, double eta,
                                              int n_zero_crossings,
                                              std::mt19937 &gen) {
  std::vector<double> x(N);
  std::uniform_real_distribution<double> dist(-eta, eta);

  do {
    for (int i = 0; i < N; ++i) {
      x[i] = dist(gen);
    }
  } while (count_zero_crossings(x) != n_zero_crossings);

  return x;
}

// ---------- Fig8 sector averages ----------
double delta_S_alpha_mean(int n_inst, double alpha, int n_therm, int n_meas,
                          int stride, double dx_width, std::mt19937 &gen) {
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;

  std::vector<double> path;
  int pinned_index = -1;
  bool enforce_orientation = false;

  if (n_inst == 0) {
    path = init_vacuum_path(N, eta, -1);
  } else if (n_inst == 1) {
    pinned_index = N / 2;
    enforce_orientation = true;
    path = init_instanton_path(N, a, eta, pinned_index);
  } else {
    return 0.0;
  }

  metropolis_alpha_constrained(path, a, eta, alpha, n_therm, dx_width, gen,
                               pinned_index, enforce_orientation);

  double sum = 0.0;
  int count = 0;

  for (int sweep = 0; sweep < n_meas; ++sweep) {
    metropolis_alpha_constrained(path, a, eta, alpha, 1, dx_width, gen,
                                 pinned_index, enforce_orientation);

    if (sweep % stride == 0) {
      const double deltaS = S_full(path, a, eta) - S_gaussian(path, a, eta);
      sum += deltaS;
      ++count;
    }
  }

  return (count > 0) ? sum / static_cast<double>(count) : 0.0;
}

double delta_S_alpha_mean_aoverride(int n_inst, double alpha, int n_therm,
                                    int n_meas, int stride, double dx_width,
                                    double a_override, std::mt19937 &gen) {
  const int N = params::N;
  const double eta = params::eta;
  const double a = a_override;

  std::vector<double> path;
  int pinned_index = -1;
  bool enforce_orientation = false;

  if (n_inst == 0) {
    path = init_vacuum_path(N, eta, -1);
  } else if (n_inst == 1) {
    pinned_index = N / 2;
    enforce_orientation = true;
    path = init_instanton_path(N, a, eta, pinned_index);
  } else {
    return 0.0;
  }

  metropolis_alpha_constrained(path, a, eta, alpha, n_therm, dx_width, gen,
                               pinned_index, enforce_orientation);

  double sum = 0.0;
  int count = 0;

  for (int sweep = 0; sweep < n_meas; ++sweep) {
    metropolis_alpha_constrained(path, a, eta, alpha, 1, dx_width, gen,
                                 pinned_index, enforce_orientation);

    if (sweep % stride == 0) {
      const double deltaS = S_full(path, a, eta) - S_gaussian(path, a, eta);
      sum += deltaS;
      ++count;
    }
  }

  return (count > 0) ? sum / static_cast<double>(count) : 0.0;
}

// ---------- Fig8 integrators ----------
double simpson_integral_aoverride(int n_alpha, int sweeps, double dx,
                                  double a_override, std::mt19937 &gen,
                                  std::ofstream *log = nullptr) {
  if (n_alpha % 2 == 0) {
    ++n_alpha;
  }

  const double h = 1.0 / (n_alpha - 1);
  double result = 0.0;

  const int n_therm = sweeps;
  const int n_meas = 20 * sweeps;
  const int stride = 10;

  for (int i = 0; i < n_alpha; ++i) {
    const double alpha =
        static_cast<double>(i) / static_cast<double>(n_alpha - 1);

    const double ds1 = delta_S_alpha_mean_aoverride(
        1, alpha, n_therm, n_meas, stride, dx, a_override, gen);
    const double ds0 = delta_S_alpha_mean_aoverride(
        0, alpha, n_therm, n_meas, stride, dx, a_override, gen);

    const double diff = ds1 - ds0;

    if (log) {
      *log << alpha << "," << ds1 << "," << ds0 << "," << diff << "\n";
    }

    if (i == 0 || i == n_alpha - 1) {
      result += diff;
    } else if (i % 2 == 1) {
      result += 4.0 * diff;
    } else {
      result += 2.0 * diff;
    }
  }

  return result * h / 3.0;
}

// ---------- misc helpers ----------
void dump_path_csv(std::ofstream &out, const std::vector<double> &x, int sector,
                   double alpha, int sample_id, double a) {
  const int N = static_cast<int>(x.size());

  for (int i = 0; i < N; ++i) {
    const double tau = i * a;
    out << sector << "," << alpha << "," << sample_id << "," << i << "," << tau
        << "," << x[i] << "\n";
  }
}

} // end anonymous namespace

double simpson_integral(int n_alpha, int sweeps, double dx, std::mt19937 &gen,
                        std::ofstream *log = nullptr) {
  if (n_alpha % 2 == 0) {
    ++n_alpha;
  }

  const double h = 1.0 / (n_alpha - 1);
  double result = 0.0;

  const int n_therm = 100;
  const int n_meas = 4000;
  const int stride = 10;

  for (int i = 0; i < n_alpha; ++i) {
    const double alpha =
        static_cast<double>(i) / static_cast<double>(n_alpha - 1);

    const double ds1 =
        delta_S_alpha_mean(1, alpha, n_therm, n_meas, stride, dx, gen);
    const double ds0 =
        delta_S_alpha_mean(0, alpha, n_therm, n_meas, stride, dx, gen);

    const double diff = ds1 - ds0;

    if (log) {
      *log << alpha << "," << ds1 << "," << ds0 << "," << diff << "\n";
    }

    if (i == 0 || i == n_alpha - 1) {
      result += diff;
    } else if (i % 2 == 1) {
      result += 4.0 * diff;
    } else {
      result += 2.0 * diff;
    }
  }

  return result * h / 3.0;
}

void run_qmidens_analysis(std::mt19937 &gen) {
  std::cout << "[📊] Running Fig.5 adiabatic switching (free energy vs T)...\n";

  const int sweeps = 1000;
  const double dx = params::dx_width;

  const int n_alpha_fine = 31;
  const int n_alpha_coarse = 15;

  const double a = params::a;

  // beta values, hence T = 1/beta
  const std::vector<double> beta_list = {40, 20, 10,  8, 6,   5,
                                         4,  3,  2.5, 2, 1.5, 1};

  // number of independent replicas per beta
  const int Nrep = 3;

  std::ofstream out("data/fig5_free_energy.csv");
  out << "beta,T,N_beta,F_harm,F_mean,F_err,I_fine_mean,I_coarse_mean,delta_"
         "logZ_mean\n";

  for (double beta : beta_list) {
    const int N_beta = static_cast<int>(std::round(beta / a));
    if (N_beta < 4)
      continue;

    std::vector<double> F_vals;
    std::vector<double> I_fine_vals;
    std::vector<double> I_coarse_vals;
    std::vector<double> delta_vals;

    F_vals.reserve(Nrep);
    I_fine_vals.reserve(Nrep);
    I_coarse_vals.reserve(Nrep);
    delta_vals.reserve(Nrep);

    std::ofstream *logptr = nullptr;
    std::ofstream log;

    for (int r = 0; r < Nrep; ++r) {
      // independent RNG stream
      std::mt19937 gen_r(gen());

      // dump integrand only for first replica at the largest beta
      if (r == 0 && std::abs(beta - beta_list.front()) < 1e-12) {
        log.open("data/fig5_switching_integrand.csv");
        log << "alpha,mean_deltaS\n";
        logptr = &log;
      } else {
        logptr = nullptr;
      }

      const double I_fine = switching_integral_free_energy_beta(
          n_alpha_fine, sweeps, dx, N_beta, gen_r, logptr);

      const double I_coarse = switching_integral_free_energy_beta(
          n_alpha_coarse, sweeps, dx, N_beta, gen_r);

      const double delta_logZ = (16.0 * I_fine - I_coarse) / 15.0;
      const double Fh = F_harmonic(beta, params::eta);
      const double Fmc = Fh + delta_logZ / beta;

      I_fine_vals.push_back(I_fine);
      I_coarse_vals.push_back(I_coarse);
      delta_vals.push_back(delta_logZ);
      F_vals.push_back(Fmc);
    }

    const double Fh = F_harmonic(beta, params::eta);

    const double F_mean = mean(F_vals);
    const double F_err = standard_error(F_vals);

    const double I_fine_mean = mean(I_fine_vals);
    const double I_coarse_mean = mean(I_coarse_vals);
    const double delta_mean = mean(delta_vals);

    const double T = 1.0 / beta;

    out << beta << "," << T << "," << N_beta << "," << Fh << "," << F_mean
        << "," << F_err << "," << I_fine_mean << "," << I_coarse_mean << ","
        << delta_mean << "\n";

    std::cout << "[beta=" << beta << "] "
              << "N_beta=" << N_beta << "  F_mean=" << F_mean
              << "  F_err=" << F_err << "\n";
  }

  std::cout << "[✓] Wrote data/fig5_free_energy.csv\n";
}

double compute_qmidens_corrected_density(int sweeps, double dx_width,
                                         int n_alpha_fine, int n_alpha_coarse,
                                         std::mt19937 &gen) {
  if (n_alpha_fine % 2 == 0) {
    ++n_alpha_fine;
  }
  if (n_alpha_coarse % 2 == 0) {
    ++n_alpha_coarse;
  }

  const double I_fine = simpson_integral(n_alpha_fine, sweeps, dx_width, gen);
  const double I_coarse =
      simpson_integral(n_alpha_coarse, sweeps, dx_width, gen);
  const double deltaS_richardson = (16.0 * I_fine - I_coarse) / 15.0;

  const double S0 = 4.0 * std::pow(params::eta, 3) / 3.0;
  const double pref = 8.0 * std::pow(params::eta, 2.5) * std::sqrt(2.0 / M_PI);

  const double n_gauss = pref * std::exp(-S0);
  const double n_corrected = n_gauss * std::exp(-deltaS_richardson);
  return n_corrected;
}

void run_fig9_paths(std::mt19937 &gen) {
  const int N = params::N;
  const double a = params::a;
  const double eta = params::eta;

  const int i0 = N / 2;

  const int n_therm = 500;
  const int stride_save = 200;
  const int n_samples = 3;

  std::vector<double> alphas = {0.0, 0.5, 1.0};

  std::ofstream out("data/fig9_paths.csv");
  out << "sector,alpha,sample_id,i,tau,x\n";

  for (double alpha : alphas) {

    // -------------------
    // 0-instanton sector
    // -------------------

    std::vector<double> x0(N, -eta);

    metropolis_alpha(x0, a, eta, alpha, n_therm, params::dx_width, gen);

    for (int s = 0; s < n_samples; ++s) {
      metropolis_alpha(x0, a, eta, alpha, stride_save, params::dx_width, gen);

      for (int i = 0; i < N; ++i) {
        double tau = i * a;
        out << 0 << "," << alpha << "," << s << "," << i << "," << tau << ","
            << x0[i] << "\n";
      }
    }

    // -------------------
    // 1-instanton sector
    // -------------------

    std::vector<double> x1(N);

    const double omega = 4.0 * eta;

    for (int i = 0; i < N; ++i) {
      double tau = (i - i0) * a;
      x1[i] = eta * std::tanh(0.5 * omega * tau);
    }

    x1[i0] = 0.0;

    metropolis_alpha(x1, a, eta, alpha, n_therm, params::dx_width, gen);

    for (int s = 0; s < n_samples; ++s) {
      metropolis_alpha(x1, a, eta, alpha, stride_save, params::dx_width, gen);

      x1[i0] = 0.0;

      for (int i = 0; i < N; ++i) {
        double tau = i * a;
        out << 1 << "," << alpha << "," << s << "," << i << "," << tau << ","
            << x1[i] << "\n";
      }
    }
  }

  std::cout << "[✓] Wrote data/fig9_paths.csv\n";
}