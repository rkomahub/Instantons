#include "core/lattice.hpp"

#include <random>

// Initialize a Euclidean lattice configuration.
Lattice::Lattice(int N, double eta, bool hot_start, std::mt19937 &gen) : x(N) {
  std::uniform_real_distribution<double> dist(-eta, eta);

  for (int i = 0; i < N; ++i) {

    // Hot start: random path. Cold start: constant vacuum configuration.
    x[i] = hot_start ? dist(gen) : eta;
  }
}

// Mutable access to one lattice site.
double &Lattice::operator[](int i) { return x[i]; }

// Read-only access to one lattice site.
const double &Lattice::operator[](int i) const { return x[i]; }

// Number of Euclidean time slices.
int Lattice::size() const { return x.size(); }

// Return the full Euclidean path.
const std::vector<double> &Lattice::get_path() const { return x; }