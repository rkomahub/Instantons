#pragma once

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
void run_basic_analysis();

/**
 * @brief Study instanton density under cooling.
 *
 * Reproduces density vs cooling sweeps,
 * revealing separation between quantum noise
 * and semiclassical instanton structure.
 */
void run_cooling_evolution_analysis();

/**
 * @brief Perform ensemble averaging.
 *
 * Computes mean correlators and statistical errors
 * for both:
 *   - Uncooled configurations
 *   - Cooled configurations
 */
void run_ensemble_analysis();

/**
 * @brief Random Instanton Liquid Model (RILM) analysis.
 *
 * Generates non-interacting multi-instanton configurations
 * and compares correlators with full Monte Carlo results.
 */
void run_rilm_analysis();

/**
 * @brief RILM with added Gaussian fluctuations (heating).
 *
 * Studies the effect of short-distance quantum fluctuations
 * on semiclassical instanton ensembles.
 */
void run_heated_rilm_analysis();

/**
 * @brief Interacting Instanton Liquid Model (IILM) analysis.
 *
 * Includes a short-range repulsive core between instantons,
 * improving semiclassical realism.
 */
void run_iilm_analysis();

/**
 * @brief Non-Gaussian correction to instanton density.
 *
 * Uses adiabatic switching to compute corrections beyond
 * the Gaussian (one-loop) approximation.
 */
void run_qmidens_analysis();