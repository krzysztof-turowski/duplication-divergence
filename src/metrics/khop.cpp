#include "khop.h"

std::vector<int> khop_degree(const Graph &G, const int &k) {
  auto shortest_paths = floyd_warshall(G);

  auto result = std::vector<int>(G.size(), 0);
  for (size_t v = 0; v < G.size(); v++) {
    for (size_t i = 0; i < G.size(); i++) {
      if (v != i && shortest_paths[v][i].distance <= k) {
        result[v]++;
      }
    }
  }

  return result;
}

std::vector<double> khop_reachability(const Graph &G, const int &k) {
  auto khop_degrees = khop_degree(G, k);

  auto result = std::vector<double>(G.size(), 0.0);
  auto degree_size = std::vector<int>(G.size(), 0);
  size_t max_degree = 0;

  for (size_t v = 0; v < G.size(); v++) {
    max_degree = std::max(max_degree, G[v].size());

    result[G[v].size()] += khop_degrees[v];
    degree_size[G[v].size()]++;
  }

  result.resize(max_degree + 1);
  for (size_t i = 0; i < result.size(); i++) {
    if (degree_size[i] > 0) {
      result[i] /= degree_size[i];
    }
  }

  return result;
}
