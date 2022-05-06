#pragma once

#include "common.h"

std::vector<double> betweenness_centrality(const SimpleGraph &);

std::vector<double> betweenness_centrality(const SimpleGraph &, const FWResults &);
