#pragma once

#include "common.h"

#include <algorithm>
#include <array>

extern const std::array<const Graph, 29> GRAPHLETS;

size_t count_graphlets(const Graph &graph, const Graph &graphlet);

size_t count_open_triangles(const Graph &G);
size_t count_triangles(const Graph &G);

size_t count_graphlets_naive(const Graph &graph, const Graph &graphlet);
