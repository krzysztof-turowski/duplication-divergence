#pragma once

#include "dd_header.h"
#include "generators/barabasi_albert.h"
#include "generators/copy_graph.h"
#include "generators/sticky.h"
#include "generators/two_step.h"
#include <random>

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

void generate_pure_dd_graph(Graph &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  for (int i = G.size(); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    G.resize(i + 1);

    for (auto j : G[parent]) {
      if (edge_distribution(generator) <= params.p) {
        G[i].insert(j), G[j].insert(i);
      }
    }
  }
}

void generate_pure_dd_connected_graph(Graph &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  for (int i = G.size(); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    G.resize(i + 1);

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
  }
}

void generate_chung_lu_graph(Graph &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  for (int i = G.size(); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    G.resize(i + 1);

    for (auto j : G[parent]) {
      if (edge_distribution(generator) <= params.p) {
        G[i].insert(j), G[j].insert(i);
      }
    }
    if (edge_distribution(generator) <= params.q) {
      G[i].insert(parent), G[parent].insert(i);
    }
  }
}

void generate_pastor_satorras_graph(Graph &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  for (int i = G.size(); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    G.resize(i + 1);

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
  }
}

void generate_graph_simple(Graph &G, const int &n, const Parameters &params) {
  switch (params.mode) {
    case Mode::STICKY:
      G = generate_sticky_graph(static_cast<const StickyParameters &>(params));
      break;

    case Mode::BA:
      generate_ba_graph(G, n, static_cast<const BarabasiAlbertParameters &>(params));
      break;

    case Mode::PURE_DUPLICATION:
      generate_pure_dd_graph(G, n, params);
      break;

    case Mode::PURE_DUPLICATION_CONNECTED:
      generate_pure_dd_connected_graph(G, n, params);
      break;

    case Mode::CHUNG_LU:
      generate_chung_lu_graph(G, n, params);
      break;

    case Mode::PASTOR_SATORRAS:
      generate_pastor_satorras_graph(G, n, params);
      break;

    case Mode::COPY_GRAPH:
      generate_copy_graph(G, n, static_cast<const CopyGraphParameters &>(params));
      break;

    case Mode::TWO_STEP:
      G = generate_two_step_graph(n, static_cast<const TwoStepParameters &>(params));
      break;

    default:
      throw std::invalid_argument("Invalid mode: " + LONG_NAME.find(params.mode)->second);
  }
}
