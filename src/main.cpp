#include "analysis_driver.hpp"

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

  run_basic_analysis();
  run_cooling_evolution_analysis();
  run_ensemble_analysis();
  run_rilm_analysis();
  run_heated_rilm_analysis();
  run_iilm_analysis();
  run_qmidens_analysis();

  return 0;
}