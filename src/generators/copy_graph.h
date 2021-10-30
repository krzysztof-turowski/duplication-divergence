#include "../dd_header.h"

class CopyGraphParameters : public Parameters {
 public:
  int a;
  int b;
  double c;

  CopyGraphParameters(const int &a, const int &b, const double &c) : a(a), b(b), c(c) { }
};

void generate_ba_graph(Graph &G, const int &n, const CopyGraphParameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  const size_t seed_size = G.size();
  G.resize(n);

  for (int i = seed_size; i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);

    const auto &a_param = std::min(params.a, i);
    while (G[i].size() < static_cast<size_t>(a_param)) {
      const int new_edge = parent_distribution(generator);
      if (G[i].find(new_edge) != G[i].end()) {
        continue;
      }
      G[i].insert(new_edge);
      G[new_edge].insert(i);
    }

    const auto &b_param = std::min(params.b, i - 1);
    std::set<int> copy_sources;
    while (copy_sources.size() < b_param) {
      const int new_copy_source = parent_distribution(generator);
      if (copy_sources.find(new_copy_source) != copy_sources.end()) {
        continue;
      }
      copy_sources.insert(new_copy_source);

      const auto &c_param = params.c;
      for (auto &&new_edge : G[new_copy_source]) {
        if (edge_distribution(generator) <= c_param && G[i].find(new_edge) == G[i].end()) {
          G[i].insert(new_edge);
          G[new_edge].insert(i);
        }
      }
    }
  }
}