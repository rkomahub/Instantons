#pragma once
#include <random>
#include <vector>

/**
 * @file analysis_driver.hpp
 * @brief High-level Monte Carlo and instanton analyses.
 *
 * These routines generate numerical observables and datasets
 * for the Euclidean double well system, including cooling,
 * correlators, and semiclassical instanton studies.
 */

/**
 * @brief Run the baseline Euclidean Monte Carlo analysis.
 *
 * Generates representative paths, correlators,
 * and instanton observables.
 */
void run_basic_analysis(std::mt19937 &gen);

/**
 * @brief Study cooling evolution of instanton configurations.
 *
 * Measures instanton density as a function
 * of the number of cooling sweeps.
 */
void run_cooling_evolution_analysis(std::mt19937 &gen);

/**
 * @brief Scan instanton density as a function of eta.
 *
 * @param etas List of eta values.
 * @param gen  Random number generator.
 */
void run_cooling_eta_scan(const std::vector<double> &etas, std::mt19937 &gen);

/**
 * @brief Compute ensemble-averaged correlators and errors.
 *
 * Includes both quantum and cooled ensembles.
 */
void run_ensemble_analysis(std::mt19937 &gen);

/**
 * @brief Run the Interacting Instanton Liquid Model analysis.
 *
 * Samples interacting instanton configurations
 * with short-range repulsion.
 */
void run_iilm_analysis(std::mt19937 &gen);