#include "io.hpp"
#include <fstream>

void save_path_to_csv(const std::vector<double> &path,
                      const std::string &filename, double a) {
  std::ofstream out(filename);
  out << "tau,x\n";
  for (size_t i = 0; i < path.size(); ++i) {
    out << i * a << "," << path[i] << "\n";
  }
}
