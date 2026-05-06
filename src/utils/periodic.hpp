#pragma once
#include <algorithm>
#include <cmath>

/**
 * @file periodic.hpp
 * @brief Small utilities for periodic one-dimensional lattices.
 */

inline int periodic_next(int i, int N) { return (i + 1 == N) ? 0 : i + 1; }

inline int periodic_prev(int i, int N) { return (i == 0) ? N - 1 : i - 1; }

// used inside iilm.cpp
inline double periodic_wrap(double t, double beta) {
  t = std::fmod(t, beta);
  if (t < 0.0) {
    t += beta;
  }
  return t;
}

inline double periodic_distance(double t1, double t2, double beta) {
  const double dt = std::fabs(t1 - t2);
  return std::min(dt, beta - dt);
}

inline double nearest_image_dt(double dt, double beta) {
  while (dt <= -0.5 * beta) {
    dt += beta;
  }

  while (dt > 0.5 * beta) {
    dt -= beta;
  }

  return dt;
}