#pragma once

#include "../graph/graph.h"
#include <cstdlib>
#include <limits>

struct FWResult {
  int32_t distance;
  int64_t number_of_paths;
};

using FWResults = std::vector<std::vector<FWResult>>;

FWResults floyd_warshall(const Graph &G);
