#include "lattice.hpp"
#include "parameters.hpp"

#include <random>

Lattice::Lattice(int N, double eta, bool hot_start) : x(N) {
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<double> dist(-eta, eta);
  for (int i = 0; i < N; ++i) {
    x[i] = hot_start ? dist(gen) : eta;
  }
}

double &Lattice::operator[](int i) { return x[i]; }
const double &Lattice::operator[](int i) const { return x[i]; }
int Lattice::size() const { return x.size(); }
const std::vector<double> &Lattice::get_path() const { return x; }
