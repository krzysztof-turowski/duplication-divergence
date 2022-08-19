#pragma once

extern "C" {
#include <nausparse.h>
#include <nauty.h>
#include <traces.h>
}

#include "common.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stack>
#include <vector>

double log_automorphisms_dense(const SimpleGraph &G);

double log_automorphisms_sparse(const SimpleGraph &G);

double log_automorphisms_traces(const SimpleGraph &G);

std::vector<int> connectivity(const SimpleGraph &G);

int isolated_nodes(const SimpleGraph &G);

int nodes_with_degree(const SimpleGraph &G, const int &d);
