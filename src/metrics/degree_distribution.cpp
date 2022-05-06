#include "degree_distribution.h"

std::vector<size_t> get_degree_distribution(const SimpleGraph &G) {
  std::vector<size_t> result(G.size(), 0);
  size_t max_degree = 0;
  for (auto &&v : G) {
    max_degree = std::max(max_degree, v.size());
    result[v.size()]++;
  }

  result.resize(max_degree + 1);
  result.shrink_to_fit();

  return result;
}
