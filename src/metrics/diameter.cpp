#include "diameter.h"

int64_t get_diameter(const Graph &G) {
  return get_diameter(G, repeated_bfs(G));
}

int64_t get_diameter(const Graph &G, const FWResults &shortest_paths) {
  auto result = 0;

  for (size_t i = 0; i < G.size(); i++) {
    for (size_t j = i + 1; j < G.size(); j++) {
      result = std::max(result, shortest_paths[i][j].distance);
    }
  }

  return result;
}
