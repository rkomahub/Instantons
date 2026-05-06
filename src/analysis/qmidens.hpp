#pragma once
#include <random>

/**
 * @file qmidens.hpp
 * @brief Non-Gaussian instanton density correction via adiabatic switching.
 *
 * The method interpolates between a Gaussian reference action
 * and the full double well action to compute free-energy corrections.
 */

/**
 * @brief Run the adiabatic switching analysis.
 *
 * Computes the switching integral using Simpson integration
 * and Richardson extrapolation, then estimates the corrected
 * instanton density.
 *
 * Intermediate integrand data are written to CSV files.
 *
 * @param gen Random number generator.
 */
void run_qmidens_analysis(std::mt19937 &gen);

/**
 * @brief Compute the non-Gaussian corrected instanton density.
 *
 * The switching integral is evaluated using Simpson integration
 * on two grids combined with Richardson extrapolation.
 *
 * @param sweeps         Metropolis sweeps per α value.
 * @param dx_width       Proposal width for Metropolis updates.
 * @param n_alpha_fine   Number of α points for the fine grid.
 * @param n_alpha_coarse Number of α points for the coarse grid.
 * @param gen            Random number generator.
 *
 * @return Corrected instanton density.
 */
double compute_qmidens_corrected_density(int sweeps, double dx_width,
                                         int n_alpha_fine, int n_alpha_coarse,
                                         std::mt19937 &gen);

/**
 * @brief Generate Monte Carlo paths used in Fig.9 (switching paths).
 *
 * Produces representative configurations in:
 *   - 0-instanton sector
 *   - 1-instanton sector (with fixed instanton location)
 *
 * The simulation samples the switched action
 *
 *     S_alpha = (1 - alpha) S_gauss + alpha S_full
 *
 * and exports configurations for plotting.
 */
void run_fig9_paths(std::mt19937 &gen);