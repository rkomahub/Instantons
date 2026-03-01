#pragma once
#include <string>

/**
 * @file ensemble.hpp
 * @brief Ensemble averaging and statistical error estimation.
 *
 * Monte Carlo observables must be averaged over independent
 * configurations in order to estimate expectation values and
 * statistical uncertainties.
 *
 * This module performs repeated simulations and computes:
 *   - Mean correlator
 *   - Standard error of the mean
 *   - Average instanton density
 */

/**
 * @brief Perform ensemble averaging over multiple Monte Carlo runs.
 *
 * For each trial:
 *   1. Generate a configuration using Metropolis sampling.
 *   2. Optionally apply cooling.
 *   3. Compute correlators.
 *   4. Count zero crossings.
 *
 * The function then computes:
 *
 *     mean( C(τ) )
 *     standard error of C(τ)
 *
 * and writes results to CSV files.
 *
 * @param trials        Number of independent configurations.
 * @param cooled        If true, apply cooling before measurement.
 * @param output_prefix Prefix used for output file names.
 */
void run_ensemble_average(int trials, bool cooled,
                          const std::string &output_prefix);