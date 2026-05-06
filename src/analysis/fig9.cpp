#include "analysis/fig9.hpp"

#include "analysis/qmidens.hpp"

#include <iostream>
#include <random>

void run_fig9_analysis(std::mt19937 &gen) {
  std::cout << "[📊] Running Fig.9 switching path generation...\n";
  run_fig9_paths(gen);
}