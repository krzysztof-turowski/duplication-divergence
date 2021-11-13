#include "betweenness.h"

struct FWResult {
  int32_t distance;
  int64_t number_of_paths;
};

std::vector<std::vector<FWResult>> floyd_warshall(const Graph &G) {
  auto result = std::vector<std::vector<FWResult>>(
      G.size(), std::vector<FWResult>(G.size(), { std::numeric_limits<int32_t>::max() / 2, 0 }));

  for (size_t i = 0; i < G.size(); i++) {
    result[i][i] = { 0, 1 };
    for (auto &&j : G[i]) {
      result[i][j] = { 1, 1 };
    }
  }

  for (size_t k = 0; k < G.size(); k++) {
    for (size_t i = 0; i < G.size(); i++) {
      for (size_t j = 0; j < G.size(); j++) {
        if (i == k || j == k) {
          continue;
        }
        if (result[i][j].distance == result[i][k].distance + result[k][j].distance) {
          result[i][j].number_of_paths +=
              result[i][k].number_of_paths * result[k][j].number_of_paths;
        }
        if (result[i][j].distance > result[i][k].distance + result[k][j].distance) {
          result[i][j].distance = result[i][k].distance + result[k][j].distance;
          result[i][j].number_of_paths =
              result[i][k].number_of_paths * result[k][j].number_of_paths;
        }
      }
    }
  }

  return result;
}

std::vector<double> betweenness_centrality(const Graph &G) {
  auto shortest_paths = floyd_warshall(G);

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
