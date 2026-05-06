#pragma once

/**
 * @file potential.hpp
 * @brief Double well potential of the Euclidean quantum system.
 *
 * The model is defined by
 *
 * \( V(x) = (x^2 - \eta^2)^2 \),
 *
 * with classical minima at \( x = \pm \eta \).
 */

/**
 * @brief Evaluate the double well potential.
 *
 * @param x    Field value.
 * @param eta  Position of the classical minima.
 *
 * @return Potential value \( V(x) \).
 */
double potential(double x, double eta);