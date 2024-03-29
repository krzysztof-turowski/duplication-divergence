#include "closeness.h"

std::vector<double> closeness(const SimpleGraph &G) {
  return closeness(G, floyd_warshall(G));
}

std::vector<double> closeness(const SimpleGraph &G, const FWResults &shortest_paths) {
  auto result = std::vector<double>(G.size(), 0.0);
  for (size_t v = 0; v < G.size(); v++) {
    int64_t distance_sum = 0;
    for (size_t i = 0; i < G.size(); i++) {
      if (shortest_paths[v][i].distance != INF) {
        distance_sum += shortest_paths[v][i].distance;
      }
    }
    result[v] = distance_sum > 0 ? static_cast<double>(G.size()) / distance_sum : 0;
  }

  return result;
}
