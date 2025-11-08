#pragma once
#include <vector>

std::vector<double> generate_iilm_path(int N, double a, double eta, int n_inst,
                                       double tau_core);