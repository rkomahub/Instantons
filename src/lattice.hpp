#pragma once
#include <random>
#include <vector>

class Lattice {
public:
  Lattice(int N, double eta, bool hot_start);

  double &operator[](int i);
  const double &operator[](int i) const;
  int size() const;
  const std::vector<double> &get_path() const;

private:
  std::vector<double> x;
};