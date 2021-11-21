#pragma once

#include "common.h"

#include <algorithm>
#include <array>

extern const std::array<const Graph, 29> GRAPHLETS;

size_t count_graphlets(const Graph &graph, const Graph &graphlet);

size_t count_open_triangles(const Graph &G);
size_t count_triangles(const Graph &G);
size_t count_four_paths(const Graph &graph);
size_t count_three_stars(const Graph &graph);
size_t count_squares(const Graph &graph);
size_t count_triangles_with_antenna(const Graph &graph);
size_t count_four_almost_cliques(const Graph &graph);
size_t count_four_cliques(const Graph &graph);
size_t count_five_paths(const Graph &graph);

size_t count_graphlets_naive(const Graph &graph, const Graph &graphlet);
