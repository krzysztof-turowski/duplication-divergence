#include "average_shortest_path.h"

double get_average_shortest_path(const Graph &G) {
  return get_average_shortest_path(G, repeated_bfs(G));
}

double get_average_shortest_path(const Graph &G, const FWResults &shortest_paths) {
  double result = 0;

  for (size_t i = 0; i < G.size(); i++) {
    for (size_t j = i + 1; j < G.size(); j++) {
      result += shortest_paths[i][j].distance;
    }
  }

  return result / (G.size() * (G.size() - 1) / 2);
}
