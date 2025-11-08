#pragma once

namespace params {
constexpr int N = 800;         // number of lattice points
constexpr double a = 0.05;     // lattice spacing
constexpr double eta = 1.4;    // potential minimum
constexpr int sweeps = 100000; // total Monte Carlo sweeps
constexpr double dx_width = 0.5;
constexpr bool hot_start = true; // true = hot, false = cold
} // namespace params