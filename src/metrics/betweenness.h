#pragma once

#include "common.h"

std::vector<double> betweenness_centrality(const Graph &);

std::vector<double> betweenness_centrality(const Graph &, const FWResults &);
