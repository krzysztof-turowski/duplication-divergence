#pragma once

#include "dd_header.h"
#include "generators/barabasi_albert.h"
#include "generators/berg.h"
#include "generators/chung_lu.h"
#include "generators/copy_graph.h"
#include "generators/kumar_linear.h"
#include "generators/pastor_satorras.h"
#include "generators/pure_dd.h"
#include "generators/pure_dd_connected.h"
#include "generators/sticky.h"
#include "generators/two_step.h"
#include <random>

SimpleGraph generate_seed_simple(const int &n0, const double &p0) {
  SimpleGraph G(n0);
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

void generate_graph_simple(SimpleGraph &G, const int &n, const Parameters &params) {
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

    case Mode::BERG:
      generate_berg_graph(G, static_cast<const BergParameters &>(params));
      break;

    case Mode::KUMAR_LINEAR:
      generate_kumar_linear_graph(G, n, static_cast<const KumarLinearParameters &>(params));
      break;

    default:
      throw std::invalid_argument("Invalid mode: " + params.long_name());
  }
}
