#pragma once

#include "../graph/graph.h"
#include <cstdlib>
#include <limits>

struct FWResult {
  int32_t distance;
  int64_t number_of_paths;
};

std::vector<std::vector<FWResult>> floyd_warshall(const Graph &G);
