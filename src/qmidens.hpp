#pragma once

void run_qmidens_analysis();

// Returns non-Gaussian corrected instanton density for current params::eta
// Uses Simpson + Richardson as in run_qmidens_analysis().
// Arguments are explicit to avoid depending on params in the header.
double compute_qmidens_corrected_density(int sweeps, double dx_width,
                                         int n_alpha_fine, int n_alpha_coarse);