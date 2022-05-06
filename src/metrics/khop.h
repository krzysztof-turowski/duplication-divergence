#pragma once

#include "common.h"

std::vector<double> khop_reachability(const SimpleGraph &G, const int &k);

std::vector<double> khop_reachability(const SimpleGraph &G, const int &k, const FWResults &);
