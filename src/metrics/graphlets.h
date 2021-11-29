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
uint64_t count_triangles_with_two_antennas(const Graph &graph);
uint64_t count_triangles_with_long_antenna(const Graph &graph);
uint64_t count_four_stars_with_edge(const Graph &graph);
uint64_t count_polygons(const Graph &graph);
uint64_t count_squares_with_antenna(const Graph &graph);
uint64_t count_almost_four_cliques_with_antenna(const Graph &graph);
uint64_t count_bowties(const Graph &graph);
uint64_t count_almost_four_cliques_with_antenna_alt(const Graph &graph);
uint64_t count_clique_two_three(const Graph &graph);
uint64_t count_houses(const Graph &graph);
uint64_t count_clique_three_one_one(const Graph &graph);
uint64_t count_four_cliques_with_antenna(const Graph &graph);
uint64_t count_three_triangles(const Graph &graph);
uint64_t count_squares_with_three_cross(const Graph &graph);
uint64_t count_four_cliques_with_flag(const Graph &graph);
uint64_t count_subdivided_crosses(const Graph &graph);
uint64_t count_five_almost_cliques(const Graph &graph);
uint64_t count_five_cliques(const Graph &graph);

uint64_t count_graphlets_naive(const Graph &graph, const Graph &graphlet);
