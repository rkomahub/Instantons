#pragma once
#include <vector>

/**
 * @file fig14_ia_interaction.hpp
 * @brief Fig.14: IA interaction in the sum ansatz + streamline via gradient
 * flow.
 */

/// Periodic minimal distance on a circle of length beta.
double periodic_distance(double t1, double t2, double beta);

/// Build sum-ansatz instanton–anti-instanton configuration on the lattice.
std::vector<double> build_sum_ansatz_ia(int N, double a, double eta,
                                        double tau_I, double tau_A);

/// Find all zero-crossings x(tau)=0 via linear interpolation.
/// Returns times in [0,beta), sorted.
std::vector<double> zero_crossings_interp(const std::vector<double> &x,
                                          double a);

/// Return tau_z = minimal periodic separation between the two crossings.
/// If crossings are not exactly 2, returns negative value.
double zero_crossing_separation(const std::vector<double> &x, double a);

/// One step of deterministic gradient flow: x <- x - eps * dS/dx.
void gradient_flow_step(std::vector<double> &x, double a, double eta,
                        double eps);