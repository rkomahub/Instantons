#pragma once
#include <random>
#include <vector>

/**
 * @file analysis_driver.hpp
 * @brief High-level analysis routines reproducing numerical experiments.
 *
 * Each function corresponds to a specific numerical experiment
 * or figure-type analysis in the instanton study of the
 * double well potential.
 *
 * These routines orchestrate:
 *   - Monte Carlo sampling
 *   - Cooling procedures
 *   - Instanton counting
 *   - Correlator computation
 *   - Instanton liquid model comparisons
 *   - Non-Gaussian density corrections
 */

/**
 * @brief Run basic Monte Carlo analysis.
 *
 * Generates:
 *   - Euclidean path configurations
 *   - Two-point correlators
 *   - Instanton counts
 *
 * Serves as the baseline full Monte Carlo study.
 */
void run_basic_analysis(std::mt19937 &gen);

/**
 * @brief Study instanton density under cooling.
 *
 * Reproduces density vs cooling sweeps,
 * revealing separation between quantum noise
 * and semiclassical instanton structure.
 */
void run_cooling_evolution_analysis(std::mt19937 &gen);

/**
 * @brief Perform a scan of the instanton density as a function of eta.
 *
 * This routine is useful for studying how the semiclassical
 * instanton density depends on the barrier height parameter.
 *
 * @param etas  List of eta values to scan.
 * @param gen   Random number generator (injected).
 */
void run_cooling_eta_scan(const std::vector<double> &etas, std::mt19937 &gen);

/**
 * @brief Perform ensemble averaging.
 *
 * Computes mean correlators and statistical errors
 * for both:
 *   - Uncooled configurations
 *   - Cooled configurations
 */
void run_ensemble_analysis(std::mt19937 &gen);

// Reproduce Schäfer Fig.7 (a,b): ensemble-averaged cooling evolution for eta
// = 1.4, 1.5, 1.6.
void run_fig7_analysis(std::mt19937 &gen);

// Reproduce Schäfer Fig.8: instanton density vs eta after cooling.
void run_fig8_analysis(const std::vector<double> &etas, std::mt19937 &gen);

void run_fig9_analysis(std::mt19937 &gen);

/**
 * @brief Random Instanton Liquid Model (RILM) analysis.
 *
 * Generates non-interacting multi-instanton configurations
 * and compares correlators with full Monte Carlo results.
 */
void run_rilm_analysis(std::mt19937 &gen);

// Reproduce Schäfer Fig.11: instanton density vs eta in the IILM.
void run_fig11_analysis(double eta);

/**
 * @brief RILM with added Gaussian fluctuations (heating).
 *
 * Studies the effect of short-distance quantum fluctuations
 * on semiclassical instanton ensembles.
 */
void run_heated_rilm_analysis(std::mt19937 &gen);

// Reproduce Schäfer Fig.14
void run_fig14_analysis(std::mt19937 &gen);

void run_fig15_analysis(std::mt19937 &gen);

// Reproduce Schäfer Fig.16: distribution of IA separations after cooling
void run_fig16_analysis(std::mt19937 &gen);

/**
 * @brief Interacting Instanton Liquid Model (IILM) analysis.
 *
 * Includes a short-range repulsive core between instantons,
 * improving semiclassical realism.
 */
void run_iilm_analysis(std::mt19937 &gen);

// Reproduce Schäfer Fig.17: instanton positions vs configuration index in IILM
void run_fig17_analysis(std::mt19937 &gen);

/**
 * @brief Non-Gaussian correction to instanton density.
 *
 * Uses adiabatic switching to compute corrections beyond
 * the Gaussian (one-loop) approximation.
 */
void run_qmidens_analysis(std::mt19937 &gen);