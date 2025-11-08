#pragma once
#include <vector>

std::vector<double> compute_correlator(const std::vector<double> &x);
void save_correlator_to_csv(const std::vector<double> &correlator, double a,
                            const std::string &filename);