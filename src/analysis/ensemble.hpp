#pragma once
#include <random>
#include <string>

/**
 * @file ensemble.hpp
 * @brief Ensemble averaging and statistical analysis.
 *
 * Computes averaged observables and statistical errors
 * from independent Monte Carlo configurations.
 */

/**
 * @brief Compute ensemble-averaged correlators and densities.
 *
 * For each configuration, the routine performs Monte Carlo sampling,
 * optional cooling, correlator measurements, and instanton counting.
 *
 * Output files contain mean values and standard errors.
 *
 * @param trials        Number of independent configurations.
 * @param cooled        Apply cooling before measurements if true.
 * @param output_prefix Prefix for output CSV files.
 * @param gen           Random number generator.
 */
void run_ensemble_average(int trials, bool cooled,
                          const std::string &output_prefix, std::mt19937 &gen);