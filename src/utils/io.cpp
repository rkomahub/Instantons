#include "utils/io.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

// Save one Euclidean path as a CSV file.
void save_path_to_csv(const std::vector<double> &path,
                      const std::string &filename, double a) {
  std::ofstream out(filename);

  out << "tau,x\n";

  // Export lattice coordinate values versus Euclidean time.
  for (size_t i = 0; i < path.size(); ++i) {
    out << i * a << "," << path[i] << "\n";
  }
}

// Save a Euclidean correlator as a CSV file.
void save_correlator_to_csv(const std::vector<double> &correlator, double a,
                            const std::string &filename) {
  std::ofstream out(filename);

  out << "tau,C(tau)\n";

  // Export correlator values at each Euclidean time separation.
  for (size_t i = 0; i < correlator.size(); ++i) {
    out << i * a << "," << correlator[i] << "\n";
  }
}

// Load a Euclidean path previously stored in CSV format.
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

    // Read only the field value x(tau).
    path.push_back(std::stod(x_str));
  }

  return path;
}