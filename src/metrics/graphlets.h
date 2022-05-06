#pragma once

#include "common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

extern const std::array<const SimpleGraph, 29> GRAPHLETS;
using AdjMatrixFn = std::function<bool(const unsigned &, const unsigned &)>;

std::array<uint64_t, GRAPHLETS.size()> count_small_graphlets(const SimpleGraph &graph);

std::array<double, GRAPHLETS.size()> relative_graphlet_frequency(
    std::array<uint64_t, GRAPHLETS.size()> const &small_graphlets);

uint64_t count_graphlets(const SimpleGraph &graph, const SimpleGraph &graphlet,
    std::function<bool(const unsigned &, const unsigned &)> *edge_exists = nullptr);

uint64_t count_open_triangles(const SimpleGraph &G, AdjMatrixFn const &adj_mat);
uint64_t count_triangles(const SimpleGraph &G, AdjMatrixFn const &adj_mat);
uint64_t count_four_paths(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_three_stars(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_squares(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_triangles_with_antenna(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_four_almost_cliques(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_four_cliques(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_five_paths(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_three_stars_with_antenna(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_four_stars(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_triangles_with_two_antennas(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_triangles_with_long_antenna(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_four_stars_with_edge(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_polygons(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_squares_with_antenna(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_almost_four_cliques_with_antenna(
    const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_bowties(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_almost_four_cliques_with_antenna_alt(
    const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_clique_two_three(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_houses(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_clique_three_one_one(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_four_cliques_with_antenna(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_three_triangles(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_squares_with_three_cross(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_four_cliques_with_flag(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_subdivided_crosses(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_five_almost_cliques(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);
uint64_t count_five_cliques(const SimpleGraph &graph, AdjMatrixFn const &adj_mat);

uint64_t count_graphlets_naive(const SimpleGraph &graph, const SimpleGraph &graphlet);
