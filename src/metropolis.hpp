#pragma once
#include "lattice.hpp"

//  Implements the lattice evolution and cooling.

class Metropolis {
public:
  Metropolis(Lattice &lattice);
  void step();

private:
  Lattice &x;
  std::mt19937 gen;
};
