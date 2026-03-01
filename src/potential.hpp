#pragma once

/**
 * @file potential.hpp
 * @brief Double well potential definition.
 *
 * The quantum mechanical model is defined by
 *
 *     V(x) = (x^2 - eta^2)^2
 *
 * The parameter eta determines:
 *   - Position of the classical minima at x = ±eta
 *   - Height of the potential barrier
 *   - Instanton action S0 = 4 eta^3 / 3
 */

/**
 * @brief Evaluate the double well potential.
 *
 * @param x    Field value.
 * @param eta  Location of classical minima.
 *
 * @return V(x) = (x^2 - eta^2)^2
 */
double potential(double x, double eta);