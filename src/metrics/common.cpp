#include "common.h"

FWResults floyd_warshall(const Graph &G) {
  auto result = FWResults(G.size(), std::vector<FWResult>(G.size(), { INF, 0 }));

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

void bfs(const Graph &graph, const unsigned &start, FWResults &results) {
  std::queue<unsigned> to_process;
  to_process.push(start);
  results[start][start].distance = 0;
  results[start][start].number_of_paths = 1;
  while (!to_process.empty()) {
    auto const &node = to_process.front();
    to_process.pop();
    for (auto &&neighbor : graph[node]) {
      if (results[start][neighbor].distance == INF) {
        results[start][neighbor].distance = results[start][node].distance + 1;
        to_process.push(neighbor);
      }
      if (results[start][neighbor].distance == results[start][node].distance + 1) {
        results[start][neighbor].number_of_paths += results[start][node].number_of_paths;
      }
    }
  }
}

FWResults repeated_bfs(const Graph &graph) {
  auto result = FWResults(graph.size(), std::vector<FWResult>(graph.size(), { INF, 0 }));
  for (size_t i = 0; i < graph.size(); i++) {
    bfs(graph, i, result);
  }
  return result;
}
