#include "graphlets.h"

using E = std::set<unsigned>;
using G = Graph;

// Corresponds to Figure 1 of "Modeling interactome: scale-free or geometric?"
// by  N. Pržulj & D. G. Corneil & I. Jurisica
const std::array<const Graph, 29> GRAPHLETS = {
  // 3-node
  G{ E{ 1 }, E{ 0, 2 }, E{ 1 } },
  G{ E{ 1, 2 }, E{ 0, 2 }, E{ 0, 1 } },
  // 4-node
  G{ E{ 1 }, E{ 0, 2 }, E{ 1, 3 }, E{ 2 } },
  G{ E{ 1, 2, 3 }, E{ 0 }, E{ 0 }, E{ 0 } },
  G{ E{ 1, 3 }, E{ 0, 2 }, E{ 1, 3 }, E{ 0, 2 } },
  G{ E{ 1, 2, 3 }, E{ 0, 2 }, E{ 0, 1 }, E{ 0 } },
  G{ E{ 1, 2, 3 }, E{ 0, 2 }, E{ 0, 1, 3 }, E{ 0, 2 } },
  G{ E{ 1, 2, 3 }, E{ 0, 2, 3 }, E{ 0, 1, 3 }, E{ 0, 1, 2 } },
  // 5-node
  G{ E{ 1 }, E{ 0, 2 }, E{ 1, 3 }, E{ 2, 4 }, E{ 3 } },
  G{ E{ 1 }, E{ 0, 2 }, E{ 1, 3, 4 }, E{ 2 }, E{ 2 } },
  G{ E{ 1, 2, 3, 4 }, E{ 0 }, E{ 0 }, E{ 0 }, E{ 0 } },
  G{ E{ 1, 2 }, E{ 0, 2, 3 }, E{ 0, 1, 4 }, E{ 1 }, E{ 2 } },
  G{ E{ 1, 2 }, E{ 0, 2 }, E{ 0, 1, 3 }, E{ 2, 4 }, E{ 3 } },
  G{ E{ 1, 2 }, E{ 0, 2 }, E{ 0, 1, 3, 4 }, E{ 2 }, E{ 2 } },
  G{ E{ 1, 4 }, E{ 0, 2 }, E{ 1, 3 }, E{ 2, 4 }, E{ 0, 3 } },
  G{ E{ 1, 3 }, E{ 0, 2 }, E{ 1, 3 }, E{ 0, 2, 4 }, E{ 3 } },
  G{ E{ 1, 2, 3 }, E{ 0, 2 }, E{ 0, 1, 3, 4 }, E{ 0, 2 }, E{ 2 } },
  G{ E{ 1, 2 }, E{ 0, 2 }, E{ 0, 1, 3, 4 }, E{ 2, 4 }, E{ 2, 3 } },
  G{ E{ 1, 2, 3 }, E{ 0, 2 }, E{ 0, 1, 3 }, E{ 0, 2, 4 }, E{ 3 } },
  G{ E{ 1, 3, 4 }, E{ 0, 2 }, E{ 1, 3, 4 }, E{ 0, 2 }, E{ 0, 2 } },
  G{ E{ 1, 3, 4 }, E{ 0, 2, 4 }, E{ 1, 3 }, E{ 0, 2 }, E{ 0, 1 } },
  G{ E{ 1, 2, 3, 4 }, E{ 0, 2 }, E{ 0, 1, 3, 4 }, E{ 0, 2 }, E{ 0, 2 } },
  G{ E{ 1, 2, 3, 4 }, E{ 0, 2, 3 }, E{ 0, 1, 3 }, E{ 0, 1, 2 }, E{ 0 } },
  G{ E{ 1, 2, 3, 4 }, E{ 0, 2, 3 }, E{ 0, 1, 4 }, E{ 0, 1 }, E{ 0, 2 } },
  G{ E{ 1, 3, 4 }, E{ 0, 2, 4 }, E{ 1, 3, 4 }, E{ 0, 2 }, E{ 0, 1, 2 } },
  G{ E{ 1, 2, 3, 4 }, E{ 0, 2, 4 }, E{ 0, 1, 3, 4 }, E{ 0, 2 }, E{ 0, 1, 2 } },
  G{ E{ 1, 3, 4 }, E{ 0, 2, 4 }, E{ 1, 3, 4 }, E{ 0, 2, 4 }, E{ 0, 1, 2, 3 } },
  G{ E{ 1, 2, 3, 4 }, E{ 0, 2, 3, 4 }, E{ 0, 1, 3, 4 }, E{ 0, 1, 2 }, E{ 0, 1, 2 } },
  G{ E{ 1, 2, 3, 4 }, E{ 0, 2, 3, 4 }, E{ 0, 1, 3, 4 }, E{ 0, 1, 2, 4 }, E{ 0, 1, 2, 3 } },
};

size_t count_graphlets(const Graph &graph, const Graph &graphlet) {
  if (graphlet == GRAPHLETS[0]) {
    return count_open_triangles(graph);
  } else if (graphlet == GRAPHLETS[1]) {
    return count_triangles(graph);
  } else if (graphlet == GRAPHLETS[2]) {
    return count_four_paths(graph);
  } else if (graphlet == GRAPHLETS[3]) {
    return count_three_stars(graph);
  } else if (graphlet == GRAPHLETS[4]) {
    return count_squares(graph);
  }
  return count_graphlets_naive(graph, graphlet);
}

size_t count_open_triangles(const Graph &graph) {
  size_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = next(it, 1); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          result++;
        }
      }
    }
  }
  return result;
}

size_t count_triangles(const Graph &graph) {
  size_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = next(it, 1); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          result++;
        }
      }
    }
  }
  return result / 3;
}

size_t count_four_paths(const Graph &graph) {
  size_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = next(it, 1); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          for (auto &&k : graph[*it]) {
            if (!v.count(k) && !graph[*jt].count(k)) {
              result++;
            }
          }
          for (auto &&k : graph[*jt]) {
            if (!v.count(k) && !graph[*it].count(k)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 2;
}

size_t count_three_stars(const Graph &graph) {
  size_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = next(it, 1); jt != v.end(); ++jt) {
        for (auto kt = next(jt, 1); kt != v.end(); ++kt) {
          if (!graph[*it].count(*jt) && !graph[*jt].count(*kt) && !graph[*kt].count(*it)) {
            result++;
          }
        }
      }
    }
  }
  return result;
}

size_t count_squares(const Graph &graph) {
  size_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = next(it, 1); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && !v.count(k) && graph[*jt].count(k)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 4;
}

bool is_isomorphic(const Graph &graph, const Graph &graphlet, std::vector<unsigned> combination) {
  do {
    for (size_t i = 0; i < combination.size() - 1; i++) {
      for (size_t j = i + 1; j < combination.size(); j++) {
        if (graphlet[i].count(j) != graph[combination[i]].count(combination[j])) {
          goto next_perm;
        }
      }
    }

    return true;
  next_perm:;
  } while (std::next_permutation(combination.begin(), combination.end()));
  return false;
}

size_t count_graphlets_naive(const Graph &graph, const Graph &graphlet) {
  std::vector<unsigned> combination(graphlet.size());
  for (size_t i = 0; i < combination.size(); i++) {
    combination[i] = i;
  }
  auto next_combination = [&]() {
    combination.back()++;
    size_t i;
    for (i = combination.size() - 1; i > 0; i--) {
      if (combination[i] >= graph.size() - (combination.size() - i - 1)) {
        combination[i - 1]++;
      } else {
        break;
      }
    }
    for (i++; i < combination.size(); i++) {
      combination[i] = combination[i - 1] + 1;
    }
    return combination.front() <= graph.size() - combination.size();
  };

  size_t result = 0;
  do {
    if (is_isomorphic(graph, graphlet, combination)) {
      result++;
    }
  } while (next_combination());
  return result;
}
