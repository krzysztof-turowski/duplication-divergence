#include "average_shortest_path.h"

double get_average_shortest_path(const SimpleGraph &G) {
  return get_average_shortest_path(G, repeated_bfs(G));
}

double get_average_shortest_path(const SimpleGraph &G, const FWResults &shortest_paths) {
  double result = 0;
  int64_t paths = 0;

  for (size_t i = 0; i < G.size(); i++) {
    for (size_t j = i + 1; j < G.size(); j++) {
      if (shortest_paths[i][j].distance != INF) {
        result += shortest_paths[i][j].distance;
        paths++;
      }
    }
  }

  return paths > 0 ? result / paths : 0;
}
