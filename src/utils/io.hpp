#pragma once
#include <string>
#include <vector>

/**
 * @file io.hpp
 * @brief CSV input/output utilities.
 */

/**
 * @brief Save a Euclidean path to CSV.
 *
 * Output columns:
 *   tau,x
 *
 * @param path      Discretized Euclidean path.
 * @param filename  Output file path.
 * @param a         Lattice spacing.
 */
void save_path_to_csv(const std::vector<double> &path,
                      const std::string &filename, double a);

/**
 * @brief Save a correlator to CSV.
 *
 * Output columns:
 *   tau,C(tau)
 *
 * @param correlator  Correlator values.
 * @param a           Lattice spacing.
 * @param filename    Output file path.
 */
void save_correlator_to_csv(const std::vector<double> &correlator, double a,
                            const std::string &filename);

/**
 * @brief Load only the path values from a CSV file with columns tau,x.
 * @param filename  Input file path.
 * @return          Vector of path values.
 */
std::vector<double> load_path_from_csv(const std::string &filename);