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

/**
 * @brief Generate Fig.7 cooling analysis data.
 */
void run_fig7_analysis(std::mt19937 &gen);

/**
 * @brief Generate Fig.8 instanton density comparison data.
 */
void run_fig8_analysis(const std::vector<double> &etas,
                       std::mt19937 &gen);

/**
 * @brief Generate Fig.9 switching-path configurations.
 */
void run_fig9_analysis(std::mt19937 &gen);

/**
 * @brief Generate Fig.10 RILM correlator data.
 */
void run_rilm_analysis(std::mt19937 &gen);

/**
 * @brief Export Gaussian effective potential data for Fig.11.
 */
void run_fig11_analysis(double eta);

/**
 * @brief Generate heated RILM data for Figs.12 and 13.
 */
void run_heated_rilm_analysis(std::mt19937 &gen);

/**
 * @brief Generate IA interaction and separation data for Figs.14 and 16.
 */
void run_fig14_analysis(std::mt19937 &gen);

/**
 * @brief Generate relaxed IA configurations for Fig.15.
 */
void run_fig15_analysis(std::mt19937 &gen);

/**
 * @brief Generate instanton trajectory data for Fig.17.
 */
void run_fig17_analysis(std::mt19937 &gen);

/**
 * @brief Run the non-Gaussian instanton density analysis.
 */
void run_qmidens_analysis(std::mt19937 &gen);