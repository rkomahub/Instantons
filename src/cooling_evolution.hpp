#pragma once
#include "lattice.hpp"
#include <string>

void run_cooling_evolution(const Lattice &original, int max_sweeps, double a,
                           const std::string &output_filename);
