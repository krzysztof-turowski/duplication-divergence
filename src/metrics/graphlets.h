#pragma once

#include "common.h"

#include <algorithm>
#include <array>
#include <cstdint>

extern const std::array<const Graph, 29> GRAPHLETS;

uint64_t count_graphlets(const Graph &graph, const Graph &graphlet);

uint64_t count_open_triangles(const Graph &G);
uint64_t count_triangles(const Graph &G);
uint64_t count_four_paths(const Graph &graph);
uint64_t count_three_stars(const Graph &graph);
uint64_t count_squares(const Graph &graph);
uint64_t count_triangles_with_antenna(const Graph &graph);
uint64_t count_four_almost_cliques(const Graph &graph);
uint64_t count_four_cliques(const Graph &graph);
uint64_t count_five_paths(const Graph &graph);
uint64_t count_three_stars_with_antenna(const Graph &graph);
uint64_t count_four_stars(const Graph &graph);

uint64_t count_graphlets_naive(const Graph &graph, const Graph &graphlet);
