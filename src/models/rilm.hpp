#pragma once
#include <random>
#include <vector>

/**
 * @file rilm.hpp
 * @brief Random Instanton Liquid Model utilities.
 *
 * Generates semiclassical configurations from randomly
 * distributed instantons and anti-instantons.
 */

/**
 * @brief Generate a random instanton liquid path.
 *
 * Builds a multi-instanton configuration from
 * independently distributed tunneling events.
 *
 * @param N        Number of lattice sites.
 * @param a        Lattice spacing.
 * @param eta      Double well parameter.
 * @param n_inst   Number of instantons and anti-instantons.
 * @param gen      Random number generator.
 *
 * @return Euclidean lattice path.
 */
std::vector<double> generate_rilm_path(int N, double a, double eta, int n_inst,
                                       std::mt19937 &gen);