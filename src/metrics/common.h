#pragma once

#include "../graph/graph.h"
#include <cstdlib>
#include <limits>
#include <queue>

auto constexpr INF = std::numeric_limits<int32_t>::max() / 2;

struct FWResult {
  int32_t distance;
  int64_t number_of_paths;
};

using FWResults = std::vector<std::vector<FWResult>>;

FWResults floyd_warshall(const SimpleGraph &G);

FWResults repeated_bfs(const SimpleGraph &G);
