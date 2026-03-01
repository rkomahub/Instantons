#pragma once
#include <random>

/**
 * @file qmidens.hpp
 * @brief Non-Gaussian correction to the instanton density via adiabatic
 * switching.
 *
 * This module implements an adiabatic switching technique to estimate
 * corrections beyond the Gaussian (one-loop) approximation to the
 * tunneling rate / instanton density.
 *
 * The idea is to interpolate between a Gaussian reference action S_gauss
 * and the full double well action S_full using:
 *
 *     S_α = (1 − α) S_gauss + α S_full ,    α ∈ [0,1]
 *
 * and compute a free-energy difference by integrating expectation values
 * along α, following the logic of Schäfer's switching method.
 */

/**
 * @brief Run the non-Gaussian correction analysis and print results.
 *
 * This function performs a full adiabatic switching computation,
 * including numerical integration (Simpson) and an error-reduced
 * estimate via Richardson extrapolation.
 *
 * It writes intermediate integrand values to CSV and prints:
 *   - I_fine, I_coarse
 *   - Richardson-extrapolated ΔS
 *   - n_gauss and n_corrected
 *
 * @param gen Random number generator (injected).
 */
void run_qmidens_analysis(std::mt19937 &gen);

/**
 * @brief Compute non-Gaussian corrected instanton density.
 *
 * The function estimates a correction factor by evaluating the integral
 * over α using Simpson's rule on two grids and Richardson extrapolation:
 *
 *     ΔS ≈ (16 I_fine − I_coarse)/15
 *
 * The corrected density is then
 *
 *     n_corrected = n_gauss * exp(−ΔS)
 *
 * where n_gauss is the semiclassical Gaussian estimate.
 *
 * @param sweeps         Metropolis sweeps used at each α point.
 * @param dx_width       Proposal width for Metropolis updates.
 * @param n_alpha_fine   Number of α points for the fine Simpson grid (odd).
 * @param n_alpha_coarse Number of α points for the coarse Simpson grid (odd).
 * @param gen            Random number generator (injected).
 *
 * @return Non-Gaussian corrected instanton density.
 */
double compute_qmidens_corrected_density(int sweeps, double dx_width,
                                         int n_alpha_fine, int n_alpha_coarse,
                                         std::mt19937 &gen);