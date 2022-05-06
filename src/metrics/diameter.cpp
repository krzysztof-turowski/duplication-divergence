#include "diameter.h"

int64_t get_diameter(const SimpleGraph &G) {
  return get_diameter(G, repeated_bfs(G));
}

int64_t get_diameter(const SimpleGraph &G, const FWResults &shortest_paths) {
  auto result = 0;

  for (size_t i = 0; i < G.size(); i++) {
    for (size_t j = i + 1; j < G.size(); j++) {
      if (shortest_paths[i][j].distance != INF) {
        result = std::max(result, shortest_paths[i][j].distance);
      }
    }
  }

  return result;
}
