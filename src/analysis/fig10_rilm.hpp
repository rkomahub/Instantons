#pragma once
#include <random>

/**
 * @brief Fig.10 Random Instanton Liquid Model correlators.
 *
 * Generates non-interacting multi-instanton configurations
 * and compares correlators with full Monte Carlo results.
 */
void run_rilm_analysis(std::mt19937 &gen);