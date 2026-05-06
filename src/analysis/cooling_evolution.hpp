#pragma once
#include "core/lattice.hpp"
#include <random>
#include <string>

/**
 * @file cooling_evolution.hpp
 * @brief Cooling analysis and instanton density measurement.
 *
 * Cooling suppresses short-distance fluctuations,
 * revealing semiclassical instanton structure.
 */

/**
 * @brief Measure instanton density during cooling evolution.
 *
 * Starting from a Monte Carlo configuration, the routine
 * performs cooling sweeps and records:
 *
 *   - instanton count
 *   - instanton density
 *   - action per instanton
 *
 * Output format:
 *
 *     n_cool, n_inst, density, action, s_per_inst
 *
 * @param original        Initial configuration.
 * @param max_sweeps      Maximum number of cooling sweeps.
 * @param a               Lattice spacing.
 * @param output_filename Output CSV filename.
 * @param gen             Random number generator.
 */
void run_cooling_evolution(const Lattice &original, int max_sweeps, double a,
                           const std::string &output_filename,
                           std::mt19937 &gen);