#pragma once
#include <string>
#include <vector>

void save_path_to_csv(const std::vector<double> &path,
                      const std::string &filename, double a);

std::vector<double> load_path_from_csv(const std::string &filename);