#include "clustering.h"

std::pair<int64_t, int64_t> count_triangles(const Graph &G) {
  std::pair<int64_t, int64_t> result{ 0, 0 };

  for (auto v : G) {
    result.second += v.size() * (v.size() - 1) / 2;

    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it, 1); jt != v.end(); ++jt) {
        if (G[*it].count(*jt)) {
          result.first++;
        }
      }
    }
  }

  result.first /= 3;

  return result;
}

double clustering_coefficient_two(const Graph &graph) {
  auto triangles = count_triangles(graph);
  return static_cast<double>(6 * triangles.first) / triangles.second;
}
