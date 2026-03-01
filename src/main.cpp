#include "analysis_driver.hpp"
#include <random>

/**
 * @file main.cpp
 * @brief Entry point of the instanton simulation project.
 *
 * This executable runs the complete numerical pipeline:
 *
 *   1. Full Monte Carlo simulation
 *   2. Cooling evolution study
 *   3. Ensemble averaging (cooled vs uncooled)
 *   4. Random Instanton Liquid Model (RILM)
 *   5. Heated RILM (Gaussian fluctuations)
 *   6. Interacting Instanton Liquid Model (IILM)
 *   7. Non-Gaussian instanton density correction
 *
 * Each routine generates CSV output files in the data/ directory.
 *
 * To run a specific analysis only, comment out unwanted calls below.
 */

int main() {

  std::mt19937 gen(std::random_device{}());
  // std::mt19937 gen(12345); // For reproducibility during development

  run_basic_analysis(gen);
  run_cooling_evolution_analysis(gen);
  run_ensemble_analysis(gen);
  run_rilm_analysis(gen);
  run_heated_rilm_analysis(gen);
  run_iilm_analysis(gen);
  run_qmidens_analysis(gen);

  return 0;
}