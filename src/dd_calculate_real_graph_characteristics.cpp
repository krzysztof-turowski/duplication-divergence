
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

void count_graph_statistics(Graph const &G, std::ostream &out) {
  const auto shortest_paths = repeated_bfs(G);

  out << "log_automorphisms_sparse: " << log_automorphisms_sparse(G) << std::endl;
  out << "get_average_shortest_path: " << get_average_shortest_path(G, shortest_paths) << std::endl;

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

  // Slower

  for (size_t i = 2; i < 5; i++) {
    const auto kh = khop_reachability(G, i, shortest_paths);
    out << "khop_reachability " << i << ": ";
    for (auto &&o : kh) {
      out << o << " ";
    }
    out << std::endl;
  }

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

int main() {
  const std::string RESULTS_FOLDER = "results/";
  const std::vector<std::string> GRAPH_NAMES = {
    "G-a-thaliana.txt",
    "G-c-elegans.txt",
    "G-college-msg.txt",
    "G-d-melanogaster.txt",
    "G-dynamic-simplewiki-10k.txt",
    "G-hep-th-citations-scc.txt",
    "G-homo-sapiens.txt",
    "G-mus-musculus.txt",
    "G-s-cerevisiae.txt",
    "G-s-pombe.txt",
  };

  std::ofstream result;
  result.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  for (auto &&graph_name : GRAPH_NAMES) {
    Graph G(read_graph_simple(FILES_FOLDER + graph_name));
    result.open(RESULTS_FOLDER + graph_name);
    count_graph_statistics(G, result);
    result.close();
  }

  return 0;
}
