#pragma once

#include "dd_header.h"

Graph generate_sticky_graph(const Parameters &params,
    std::uniform_real_distribution<double> &distribution, std::mt19937 &generator) {
  const auto &degrees = params.degrees;
  Graph G(degrees.size());
  int64_t degree_sum = 0;
  for (auto &&deg : params.degrees) {
    degree_sum += deg;
  }

  for (size_t i = 0; i < G.size(); i++) {
    for (size_t j = 0; j < G.size(); j++) {
      if (distribution(generator) * degree_sum * degree_sum
          <= static_cast<double>(degrees[i]) * degrees[j]) {
        G[i].insert(j), G[j].insert(i);
      }
    }
  }

  return G;
}

Graph generate_ba_graph(Graph &&G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());

  const uint32_t &m_param = params.m;

  const size_t seed_size = G.size();
  if (m_param > seed_size) {
    throw std::invalid_argument("m parameter must not be greater than seed graph size.");
  }
  G.resize(n);

  std::vector<int> degrees(n, 0);
  for (size_t i = 0; i < seed_size; i++) {
    degrees[i] = G[i].size();
  }

  for (int i = seed_size; i < n; i++) {
    std::discrete_distribution<int> distribution(degrees.begin(), degrees.end());
    while (G[i].size() < m_param) {
      const int new_edge = distribution(generator);
      if (G[i].find(new_edge) != G[i].end()) {
        continue;
      }

      degrees[i]++;
      degrees[new_edge]++;
      G[i].insert(new_edge);
      G[new_edge].insert(i);
    }
  }
  return std::move(G);
}

Graph generate_seed_simple(const int &n0, const double &p0) {
  Graph G(n0);
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  for (int i = 0; i < n0; i++) {
    for (int j = i + 1; j < n0; j++) {
      if (edge_distribution(generator) <= p0) {
        G[i].insert(j), G[j].insert(i);
      }
    }
  }
  return G;
}

void generate_graph_simple(Graph &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  if (params.mode == Mode::STICKY) {
    G = generate_sticky_graph(params, edge_distribution, generator);
    return;
  }
  if (params.mode == Mode::BA) {
    G = generate_ba_graph(std::move(G), n, params);
    return;
  }

  for (int i = G.size(); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    G.resize(i + 1);
    if (params.mode == Mode::PURE_DUPLICATION) {
      for (auto j : G[parent]) {
        if (edge_distribution(generator) <= params.p) {
          G[i].insert(j), G[j].insert(i);
        }
      }
    } else if (params.mode == Mode::PURE_DUPLICATION_CONNECTED) {
      while (true) {
        for (auto j : G[parent]) {
          if (edge_distribution(generator) <= params.p) {
            G[i].insert(j), G[j].insert(i);
          }
        }
        if (!G[i].empty()) {
          break;
        }
        parent = parent_distribution(generator);
      }
    } else if (params.mode == Mode::CHUNG_LU) {
      for (auto j : G[parent]) {
        if (edge_distribution(generator) <= params.p) {
          G[i].insert(j), G[j].insert(i);
        }
      }
      if (edge_distribution(generator) <= params.q) {
        G[i].insert(parent), G[parent].insert(i);
      }
    } else if (params.mode == Mode::PASTOR_SATORRAS) {
      for (int j = 0; j < i; j++) {
        if (G[parent].count(j)) {
          if (edge_distribution(generator) <= params.p) {
            G[i].insert(j), G[j].insert(i);
          }
        } else {
          if (edge_distribution(generator) <= params.r / i) {
            G[i].insert(j), G[j].insert(i);
          }
        }
      }
    } else {
      throw std::invalid_argument("Invalid mode: " + LONG_NAME.find(params.mode)->second);
    }
  }
  return;
}
