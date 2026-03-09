#include "io.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

void save_path_to_csv(const std::vector<double> &path,
                      const std::string &filename, double a) {
  std::ofstream out(filename);
  out << "tau,x\n";
  for (size_t i = 0; i < path.size(); ++i) {
    out << i * a << "," << path[i] << "\n";
  }
}

std::vector<double> load_path_from_csv(const std::string &filename) {
  std::ifstream in(filename);
  if (!in) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  std::vector<double> path;
  std::string line;

  // Skip header
  std::getline(in, line);

  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string tau_str, x_str;

    std::getline(ss, tau_str, ',');
    std::getline(ss, x_str, ',');

    path.push_back(std::stod(x_str));
  }

  return path;
}