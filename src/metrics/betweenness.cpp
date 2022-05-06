#include "betweenness.h"

std::vector<double> betweenness_centrality(const SimpleGraph &G) {
  return betweenness_centrality(G, floyd_warshall(G));
}

std::vector<double> betweenness_centrality(const SimpleGraph &G, const FWResults &shortest_paths) {
  auto result = std::vector<double>(G.size(), 0.0);
  for (size_t v = 0; v < G.size(); v++) {
    for (size_t i = 0; i < G.size(); i++) {
      for (size_t j = 0; j < G.size(); j++) {
        const auto &ij = shortest_paths[i][j];
        const auto &iv = shortest_paths[i][v];
        const auto &vj = shortest_paths[v][j];
        if (ij.number_of_paths != 0 && i != v && v != j && j != i
            && iv.distance + vj.distance == ij.distance) {
          result[v] +=
              static_cast<double>(iv.number_of_paths * vj.number_of_paths) / ij.number_of_paths;
        }
      }
    }
  }

  return result;
}
