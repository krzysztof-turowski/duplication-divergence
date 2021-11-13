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

double log_automorphisms_dense(const Graph &G);

double log_automorphisms_sparse(const Graph &G);

double log_automorphisms_traces(const Graph &G);

std::vector<int> connectivity(const Graph &G);

int isolated_nodes(const Graph &G);
