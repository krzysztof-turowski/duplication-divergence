#pragma once

#include "./dd_dag.h"

#include <gmpxx.h>

#include <algorithm>
#include <map>
#include <random>
#include <set>
#include <utility>
#include <vector>

typedef std::pair<int, int> VertexPair;

const int AGE_ZERO = 0;
const int AGE_MAX = std::numeric_limits<int>::max();

std::set<VertexPair> get_perfect_pairs(
    const std::vector<int> &node_age, const double &fraction) {
  if (std::isnan(fraction)) {
    return std::set<VertexPair>();
  }
  std::vector<VertexPair> count_age;
  for (size_t i = 0; i < node_age.size(); i++) {
    if (node_age[i] == AGE_ZERO) {
      continue;
    }
    for (size_t j = i + 1; j < node_age.size(); j++) {
      if (node_age[j] == AGE_ZERO) {
        continue;
      }
      if (node_age[i] < node_age[j]) {
        count_age.push_back(std::make_pair(i, j));
      } else if (node_age[j] < node_age[i]) {
        count_age.push_back(std::make_pair(j, i));
      }
    }
  }
  std::random_device device;
  std::mt19937 generator(device());
  std::set<VertexPair> perfect_pairs;
  for (size_t i = count_age.size() - 1; i >= (1.0 - fraction) * count_age.size(); i--) {
    std::uniform_int_distribution<int> swap_distribution(0, i);
    int index = swap_distribution(generator);
    std::swap(count_age[i], count_age[index]);
    perfect_pairs.insert(count_age[i]);
  }
  return perfect_pairs;
}

DAG get_DAG_from_perfect_pairs(const std::set<VertexPair> &perfect_pairs, const int &n) {
  DAG G(n);
  for (const auto &uv : perfect_pairs) {
    G.add_edge(uv.second, uv.first);
  }
  return G;
}

void set_perfect_pairs(
    std::map<VertexPair, long double> p_uv, const std::set<VertexPair> &perfect_pairs) {
  for (const auto &uv : perfect_pairs) {
    const auto vu(std::make_pair(uv.second, uv.first));
    p_uv.insert(std::make_pair(uv, 1.0L)), p_uv.insert(std::make_pair(vu, 0.0L));
  }
}
