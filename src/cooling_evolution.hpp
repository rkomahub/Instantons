#pragma once
#include "lattice.hpp"
#include <string>

/**
 * @file cooling_evolution.hpp
 * @brief Study of instanton density under cooling.
 *
 * Cooling is a procedure in which only updates that decrease
 * the Euclidean action are accepted.
 *
 * This removes short-distance quantum fluctuations and reveals
 * semiclassical structures such as instantons.
 */

/**
 * @brief Perform cooling sweeps and monitor instanton density.
 *
 * Starting from an original Monte Carlo configuration,
 * the function performs successive cooling sweeps and measures:
 *
 *   - Number of zero crossings N_{I+A}
 *   - Instanton density n_{I+A} = N_{I+A} / beta
 *
 * as a function of the number of cooling sweeps.
 *
 * The output CSV file contains:
 *
 *     n_cool, n_inst, density
 *
 * where:
 *   n_cool  = number of cooling sweeps performed
 *   n_inst  = number of zero crossings
 *   density = n_inst / beta
 *
 * @param original        Initial Monte Carlo configuration.
 * @param max_sweeps      Maximum number of cooling sweeps.
 * @param a               Lattice spacing.
 * @param output_filename Output CSV file path.
 */
void run_cooling_evolution(const Lattice &original, int max_sweeps, double a,
                           const std::string &output_filename);