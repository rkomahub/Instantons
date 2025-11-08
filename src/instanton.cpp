#include "instanton.hpp"

int count_zero_crossings(const std::vector<double> &path) {
  int count = 0;
  int N = path.size();

  for (int i = 0; i < N; ++i) {
    int j = (i + 1) % N;
    if ((path[i] < 0 && path[j] > 0) || (path[i] > 0 && path[j] < 0)) {
      ++count;
    }
  }
  return count;
}
