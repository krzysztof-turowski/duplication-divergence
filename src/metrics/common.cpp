#include "common.h"

FWResults floyd_warshall(const Graph &G) {
  auto result = FWResults(
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
