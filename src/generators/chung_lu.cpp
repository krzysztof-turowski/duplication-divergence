#include "chung_lu.h"

void generate_chung_lu_graph(SimpleGraph &G, const int &n, const Parameters &params) {
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
