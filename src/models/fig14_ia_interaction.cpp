#include "models/fig14_ia_interaction.hpp"
#include "core/potential.hpp"

#include <algorithm>
#include <cmath>

double periodic_distance(double t1, double t2, double beta) {
  const double dt = std::fabs(t1 - t2);
  return std::min(dt, beta - dt);
}

static double wrap_to_nearest_image(double dt, double beta) {
  // map dt to (-beta/2, beta/2]
  while (dt <= -0.5 * beta)
    dt += beta;
  while (dt > +0.5 * beta)
    dt -= beta;
  return dt;
}

std::vector<double> build_sum_ansatz_ia(int N, double a, double eta,
                                        double tau_I, double tau_A) {
  const double omega = 4.0 * eta;

  std::vector<double> x(N, 0.0);

  for (int i = 0; i < N; ++i) {
    const double tau = i * a;

    // For Fig.14 the pair is far from the boundaries, so use the
    // non-wrapped infinite-line ansatz to avoid artificial branch cuts.
    const double dI = tau - tau_I;
    const double dA = tau - tau_A;

    const double tI = std::tanh(0.5 * omega * dI);
    const double tA = std::tanh(0.5 * omega * dA);

    x[i] = eta * (tI - tA - 1.0);
  }

  return x;
}

std::vector<double> zero_crossings_interp(const std::vector<double> &x,
                                          double a) {
  const int N = static_cast<int>(x.size());
  const double beta = N * a;

  std::vector<double> zeros;
  zeros.reserve(8);

  for (int i = 0; i < N; ++i) {
    const int j = (i + 1) % N;

    const double xi = x[i];
    const double xj = x[j];

    // If exactly zero at a site, record it.
    if (xi == 0.0) {
      zeros.push_back(i * a);
      continue;
    }

    const bool sign_change = (xi < 0.0 && xj > 0.0) || (xi > 0.0 && xj < 0.0);
    if (!sign_change)
      continue;

    // Linear interpolation on the link (i -> i+1). Handle last link as [tau_i,
    // beta].
    const double tau_i = i * a;
    const double tau_j = (i == N - 1) ? beta : (i + 1) * a;

    // Solve xi + t*(xj - xi) = 0 => t = xi / (xi - xj)
    const double t = xi / (xi - xj);
    double tau0 = tau_i + t * (tau_j - tau_i);

    // Map beta -> 0
    if (tau0 >= beta)
      tau0 -= beta;
    if (tau0 < 0.0)
      tau0 += beta;

    zeros.push_back(tau0);
  }

  std::sort(zeros.begin(), zeros.end());
  // Remove near-duplicates (can happen if xi==0 and sign-change in adjacent
  // link)
  zeros.erase(std::unique(zeros.begin(), zeros.end(),
                          [](double a1, double a2) {
                            return std::fabs(a1 - a2) < 1e-12;
                          }),
              zeros.end());

  return zeros;
}

double zero_crossing_separation(const std::vector<double> &x, double a) {
  const int N = static_cast<int>(x.size());
  const double beta = N * a;

  auto z = zero_crossings_interp(x, a);
  if (z.size() != 2)
    return -1.0;

  return periodic_distance(z[0], z[1], beta);
}

void gradient_flow_step(std::vector<double> &x, double a, double eta,
                        double eps) {
  const int N = static_cast<int>(x.size());
  std::vector<double> grad(N, 0.0);

  // dS/dx_i = (2x_i - x_{i-1} - x_{i+1})/(2a) + a * dV/dx_i
  // V = (x^2 - eta^2)^2 => dV/dx = 4x(x^2 - eta^2)
  for (int i = 0; i < N; ++i) {
    const int ip = (i + 1) % N;
    const int im = (i - 1 + N) % N;

    const double lap = (2.0 * x[i] - x[im] - x[ip]) / (2.0 * a);
    const double dV = 4.0 * x[i] * (x[i] * x[i] - eta * eta);

    grad[i] = lap + a * dV;
  }

  for (int i = 0; i < N; ++i) {
    x[i] -= eps * grad[i];
  }
}