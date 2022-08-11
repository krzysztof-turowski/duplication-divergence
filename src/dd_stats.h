#include "./dd_input.h"
#include "./metrics/automorphisms.h"
#include "./metrics/average_shortest_path.h"
#include "./metrics/betweenness.h"
#include "./metrics/closeness.h"
#include "./metrics/clustering.h"
#include "./metrics/common.h"
#include "./metrics/degree_distribution.h"
#include "./metrics/diameter.h"
#include "./metrics/graphlets.h"
#include "./metrics/khop.h"

#include <algorithm>
#include <chrono>
#include <string>

constexpr int FAST = 1 << 0;
constexpr int SLOW = 1 << 1;

void count_graph_statistics(SimpleGraph const &G, std::ostream &out, int const &mode) {
  if (mode & FAST) {
    const auto shortest_paths = repeated_bfs(G);

    out << "log_automorphisms_sparse: " << log_automorphisms_sparse(G) << std::endl;
    out << "get_average_shortest_path: " << get_average_shortest_path(G, shortest_paths)
        << std::endl;

    {
      const auto bc = betweenness_centrality(G, shortest_paths);
      out << "betweenness_centrality: ";
      for (auto &&i : bc) {
        out << i << " ";
      }
      out << std::endl;
    }
    {
      const auto cl = closeness(G, shortest_paths);
      out << "closeness: ";
      for (auto &&i : cl) {
        out << i << " ";
      }
      out << std::endl;
    }

    out << "clustering_coefficient_two: " << clustering_coefficient_two(G) << std::endl;
    out << "clustering_coefficient_three: " << clustering_coefficient_three(G) << std::endl;

    {
      const auto cl = get_degree_distribution(G);
      out << "get_degree_distribution: ";
      for (auto &&i : cl) {
        out << i << " ";
      }
      out << std::endl;
    }

    out << "get_diameter: " << get_diameter(G, shortest_paths) << std::endl;
    for (size_t i = 2; i < 5; i++) {
      const auto kh = khop_reachability(G, i, shortest_paths);
      out << "khop_reachability " << i << ": ";
      for (auto &&o : kh) {
        out << o << " ";
      }
      out << std::endl;
    }
  }

  if (mode & SLOW) {
    out << "clustering_coefficient_four: " << clustering_coefficient_four(G) << std::endl;

    auto small_graphlets = count_small_graphlets(G);
    out << "Graphlets: ";
    for (auto &&graphlet : small_graphlets) {
      out << graphlet << " ";
    }
    out << std::endl;

    out << "RGF: ";
    auto rgfs = relative_graphlet_frequency(small_graphlets);
    for (auto &&rgf : rgfs) {
      out << rgf << " ";
    }
    out << std::endl;
  }
}
