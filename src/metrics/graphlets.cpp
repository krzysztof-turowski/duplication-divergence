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

uint64_t count_graphlets(const Graph &graph, const Graph &graphlet, AdjMatrixFn *edge_exists) {
  auto const &func =
      edge_exists == nullptr ? [&](auto u, auto v) -> bool { return graph[u].count(v); }
                             : *edge_exists;

  if (graphlet == GRAPHLETS[0]) {
    return count_open_triangles(graph, func);
  } else if (graphlet == GRAPHLETS[1]) {
    return count_triangles(graph, func);
  } else if (graphlet == GRAPHLETS[2]) {
    return count_four_paths(graph, func);
  } else if (graphlet == GRAPHLETS[3]) {
    return count_three_stars(graph, func);
  } else if (graphlet == GRAPHLETS[4]) {
    return count_squares(graph, func);
  } else if (graphlet == GRAPHLETS[5]) {
    return count_triangles_with_antenna(graph, func);
  } else if (graphlet == GRAPHLETS[6]) {
    return count_four_almost_cliques(graph, func);
  } else if (graphlet == GRAPHLETS[7]) {
    return count_four_cliques(graph, func);
  } else if (graphlet == GRAPHLETS[8]) {
    return count_five_paths(graph, func);
  } else if (graphlet == GRAPHLETS[9]) {
    return count_three_stars_with_antenna(graph, func);
  } else if (graphlet == GRAPHLETS[10]) {
    return count_four_stars(graph, func);
  } else if (graphlet == GRAPHLETS[11]) {
    return count_triangles_with_two_antennas(graph, func);
  } else if (graphlet == GRAPHLETS[12]) {
    return count_triangles_with_long_antenna(graph, func);
  } else if (graphlet == GRAPHLETS[13]) {
    return count_four_stars_with_edge(graph, func);
  } else if (graphlet == GRAPHLETS[14]) {
    return count_polygons(graph, func);
  } else if (graphlet == GRAPHLETS[15]) {
    return count_squares_with_antenna(graph, func);
  } else if (graphlet == GRAPHLETS[16]) {
    return count_almost_four_cliques_with_antenna(graph, func);
  } else if (graphlet == GRAPHLETS[17]) {
    return count_bowties(graph, func);
  } else if (graphlet == GRAPHLETS[18]) {
    return count_almost_four_cliques_with_antenna_alt(graph, func);
  } else if (graphlet == GRAPHLETS[19]) {
    return count_clique_two_three(graph, func);
  } else if (graphlet == GRAPHLETS[20]) {
    return count_houses(graph, func);
  } else if (graphlet == GRAPHLETS[21]) {
    return count_clique_three_one_one(graph, func);
  } else if (graphlet == GRAPHLETS[22]) {
    return count_four_cliques_with_antenna(graph, func);
  } else if (graphlet == GRAPHLETS[23]) {
    return count_three_triangles(graph, func);
  } else if (graphlet == GRAPHLETS[24]) {
    return count_squares_with_three_cross(graph, func);
  } else if (graphlet == GRAPHLETS[25]) {
    return count_four_cliques_with_flag(graph, func);
  } else if (graphlet == GRAPHLETS[26]) {
    return count_subdivided_crosses(graph, func);
  } else if (graphlet == GRAPHLETS[27]) {
    return count_five_almost_cliques(graph, func);
  }
  return count_graphlets_naive(graph, graphlet);
}

uint64_t count_open_triangles(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          result++;
        }
      }
    }
  }
  return result;
}

uint64_t count_triangles(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (adj_mat(*it, *jt)) {
          result++;
        }
      }
    }
  }
  return result / 3;
}

uint64_t count_four_paths(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          for (auto &&k : graph[*it]) {
            if (!adj_mat(i, k) && !adj_mat(*jt, k)) {
              result++;
            }
          }
          for (auto &&k : graph[*jt]) {
            if (!adj_mat(i, k) && !adj_mat(*it, k)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 2;
}

uint64_t count_three_stars(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!adj_mat(*it, *jt) && !adj_mat(*jt, *kt) && !adj_mat(*kt, *it)) {
            result++;
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_squares(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && !adj_mat(i, k) && adj_mat(*jt, k)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 4;
}

uint64_t count_triangles_with_antenna(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (adj_mat(*it, *jt) + adj_mat(*jt, *kt) + adj_mat(*kt, *it) == 1) {
            result++;
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_almost_cliques(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (adj_mat(*it, *jt) + adj_mat(*jt, *kt) + adj_mat(*kt, *it) == 2) {
            result++;
          }
        }
      }
    }
  }
  return result / 2;
}

uint64_t count_four_cliques(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.upper_bound(i); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (adj_mat(*it, *jt) + adj_mat(*jt, *kt) + adj_mat(*kt, *it) == 3) {
            result++;
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_five_paths(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          std::vector<unsigned> ks;
          for (auto &&k : graph[*it]) {
            if (!adj_mat(i, k) && !adj_mat(*jt, k)) {
              ks.push_back(k);
            }
          }
          for (auto &&l : graph[*jt]) {
            if (!adj_mat(i, l) && !adj_mat(*it, l)) {
              for (auto &&k : ks) {
                if (!adj_mat(l, k)) {
                  result++;
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_three_stars_with_antenna(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (adj_mat(*it, *jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!adj_mat(*jt, *kt) && !adj_mat(*kt, *it)) {
            for (auto &&l : graph[*it]) {
              if (!adj_mat(i, l) && !adj_mat(*jt, l) && !adj_mat(*kt, l)) {
                result++;
              }
            }

            for (auto &&l : graph[*jt]) {
              if (!adj_mat(i, l) && !adj_mat(*it, l) && !adj_mat(*kt, l)) {
                result++;
              }
            }

            for (auto &&l : graph[*kt]) {
              if (!adj_mat(i, l) && !adj_mat(*it, l) && !adj_mat(*jt, l)) {
                result++;
              }
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_stars(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (adj_mat(*it, *jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (adj_mat(*jt, *kt) || adj_mat(*it, *kt)) {
            continue;
          }

          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (!adj_mat(*kt, *lt) && !adj_mat(*lt, *it) && !adj_mat(*jt, *lt)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_triangles_with_two_antennas(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          continue;
        }

        for (auto &&k : graph[*it]) {
          if (adj_mat(*jt, k) != 0 || adj_mat(i, k) != 0) {
            continue;
          }

          for (auto &&l : graph[*jt]) {
            if (adj_mat(l, k) == 0 && adj_mat(l, *it) == 0 && adj_mat(l, i) == 0) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_triangles_with_long_antenna(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (adj_mat(*it, *jt)) {
          continue;
        }

        for (auto kt = graph[*it].begin(); kt != graph[*it].end(); ++kt) {
          if (adj_mat(*kt, *jt) || adj_mat(*kt, i)) {
            continue;
          }
          for (auto lt = std::next(kt); lt != graph[*it].end(); ++lt) {
            if (!adj_mat(*lt, *jt) && !adj_mat(*lt, i) && adj_mat(*lt, *kt)) {
              result++;
            }
          }
        }

        for (auto kt = graph[*jt].begin(); kt != graph[*jt].end(); ++kt) {
          if (adj_mat(*kt, *it) || adj_mat(*kt, i)) {
            continue;
          }
          for (auto lt = std::next(kt); lt != graph[*jt].end(); ++lt) {
            if (!adj_mat(*lt, *it) && !adj_mat(*lt, i) && adj_mat(*lt, *kt)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_stars_with_edge(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (adj_mat(*kt, *it) + adj_mat(*kt, *jt) + adj_mat(*kt, *lt) + adj_mat(*lt, *it)
                    + adj_mat(*lt, *jt) + adj_mat(*it, *jt)
                == 1) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_polygons(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto &&j : v) {
      for (auto &&k : graph[j]) {
        if (adj_mat(i, k)) {
          continue;
        }
        for (auto &&l : graph[k]) {
          if (adj_mat(i, l) || adj_mat(j, l)) {
            continue;
          }
          for (auto &&m : graph[l]) {
            if (adj_mat(i, m) && !adj_mat(j, m) && !adj_mat(k, m)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 10;
}

uint64_t count_squares_with_antenna(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && !adj_mat(i, k) && adj_mat(*jt, k)) {
              for (auto &&l : v) {
                if (l != *it && l != *jt && l != k && !adj_mat(l, *it) && !adj_mat(l, *jt)
                    && !adj_mat(l, k)) {
                  result++;
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_almost_four_cliques_with_antenna(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && adj_mat(i, k) && adj_mat(*jt, k)) {
              for (auto &&l : v) {
                if (l != *it && l != *jt && l != k && !adj_mat(l, *it) && !adj_mat(l, *jt)
                    && !adj_mat(l, k)) {
                  result++;
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_bowties(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (adj_mat(*it, *jt) && adj_mat(*kt, *lt) && !adj_mat(*it, *kt) && !adj_mat(*it, *lt)
                && !adj_mat(*jt, *kt) && !adj_mat(*jt, *lt)) {
              result++;
            } else if (adj_mat(*it, *kt) && adj_mat(*jt, *lt) && !adj_mat(*it, *jt)
                       && !adj_mat(*it, *lt) && !adj_mat(*kt, *jt) && !adj_mat(*kt, *lt)) {
              result++;
            } else if (adj_mat(*it, *lt) && adj_mat(*jt, *kt) && !adj_mat(*it, *jt)
                       && !adj_mat(*it, *kt) && !adj_mat(*lt, *jt) && !adj_mat(*lt, *kt)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_almost_four_cliques_with_antenna_alt(
    const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (adj_mat(*it, *jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && !adj_mat(i, k) && adj_mat(*jt, k)) {
              for (auto &&l : v) {
                if (l != *it && l != *jt && l != k && !adj_mat(l, *it) && !adj_mat(l, *jt)
                    && !adj_mat(l, k)) {
                  result++;
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_clique_two_three(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (adj_mat(*it, *jt)) {
          continue;
        }

        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (adj_mat(*it, *kt) || adj_mat(*jt, *kt)) {
            continue;
          }

          for (auto &&l : graph[*it]) {
            if (i < l && l != *jt && l != *kt && adj_mat(*jt, l) && adj_mat(*kt, l)
                && !adj_mat(i, l)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_houses(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          continue;
        }

        for (auto &&k : graph[*it]) {
          if (adj_mat(*jt, k) || adj_mat(i, k)) {
            continue;
          }

          for (auto &&l : graph[*jt]) {
            if (adj_mat(l, k) && !adj_mat(l, *it) && !adj_mat(l, i)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_clique_three_one_one(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (adj_mat(*it, *jt)) {
          continue;
        }

        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (adj_mat(*it, *kt) || adj_mat(*jt, *kt)) {
            continue;
          }

          for (auto &&l : graph[*it]) {
            if (i < l && l != *jt && l != *kt && adj_mat(*jt, l) && adj_mat(*kt, l)
                && adj_mat(i, l)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_cliques_with_antenna(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.upper_bound(i); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!adj_mat(*jt, *kt) || !adj_mat(*kt, *it)) {
            continue;
          }

          for (auto &&l : v) {
            if (l != *it && l != *jt && l != *kt && !adj_mat(*it, l) && !adj_mat(*jt, l)
                && !adj_mat(*kt, l)) {
              result++;
            }
          }
          for (auto &&l : graph[*it]) {
            if (l != i && l != *jt && l != *kt && !adj_mat(i, l) && !adj_mat(*jt, l)
                && !adj_mat(*kt, l)) {
              result++;
            }
          }
          for (auto &&l : graph[*jt]) {
            if (l != *it && l != i && l != *kt && !adj_mat(*it, l) && !adj_mat(i, l)
                && !adj_mat(*kt, l)) {
              result++;
            }
          }
          for (auto &&l : graph[*kt]) {
            if (l != *it && l != *jt && l != i && !adj_mat(*it, l) && !adj_mat(*jt, l)
                && !adj_mat(i, l)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_three_triangles(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;

  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            auto id = adj_mat(*it, *jt) + adj_mat(*it, *kt) + adj_mat(*it, *lt);
            auto jd = adj_mat(*it, *jt) + adj_mat(*jt, *kt) + adj_mat(*jt, *lt);
            auto kd = adj_mat(*kt, *jt) + adj_mat(*it, *kt) + adj_mat(*kt, *lt);
            auto ld = adj_mat(*lt, *jt) + adj_mat(*lt, *kt) + adj_mat(*it, *lt);
            if (id + jd + kd + ld == 6 && (id == 1 || jd == 1 || kd == 1 || ld == 1)
                && (id != 3 && jd != 3 && kd != 3 && ld != 3)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_squares_with_three_cross(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;

  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (i < *it && adj_mat(*it, *jt) && adj_mat(*it, *kt) && !adj_mat(*jt, *kt)) {
            for (auto &&l : graph[*jt]) {
              if (!adj_mat(i, l) && l != *it && l != *kt && l != i && !adj_mat(*it, l)
                  && adj_mat(*kt, l)) {
                result++;
              }
            }
          }

          if (i < *jt && adj_mat(*jt, *it) && adj_mat(*jt, *kt) && !adj_mat(*it, *kt)) {
            for (auto &&l : graph[*kt]) {
              if (!adj_mat(i, l) && l != *jt && l != *it && l != i && !adj_mat(*jt, l)
                  && adj_mat(*it, l)) {
                result++;
              }
            }
          }

          if (i < *kt && adj_mat(*kt, *jt) && adj_mat(*kt, *it) && !adj_mat(*jt, *it)) {
            for (auto &&l : graph[*it]) {
              if (!adj_mat(i, l) && l != *kt && l != *jt && l != i && !adj_mat(*kt, l)
                  && adj_mat(*jt, l)) {
                result++;
              }
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_cliques_with_flag(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.upper_bound(i); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!adj_mat(*jt, *kt) || !adj_mat(*kt, *it)) {
            continue;
          }

          for (auto &&l : v) {
            if (l != *it && l != *jt && l != *kt
                && (adj_mat(*it, l) + adj_mat(*jt, l) + adj_mat(*kt, l) == 1)) {
              result++;
            }
          }
          for (auto &&l : graph[*it]) {
            if (l != i && l != *jt && l != *kt
                && (adj_mat(i, l) + adj_mat(*jt, l) + adj_mat(*kt, l) == 1)) {
              result++;
            }
          }
          for (auto &&l : graph[*jt]) {
            if (l != *it && l != i && l != *kt
                && (adj_mat(*it, l) + adj_mat(i, l) + adj_mat(*kt, l) == 1)) {
              result++;
            }
          }
          for (auto &&l : graph[*kt]) {
            if (l != *it && l != *jt && l != i
                && (adj_mat(*it, l) + adj_mat(*jt, l) + adj_mat(i, l) == 1)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 2;
}

uint64_t count_subdivided_crosses(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;

  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            auto id = adj_mat(*it, *jt) + adj_mat(*it, *kt) + adj_mat(*it, *lt);
            auto jd = adj_mat(*it, *jt) + adj_mat(*jt, *kt) + adj_mat(*jt, *lt);
            auto kd = adj_mat(*kt, *jt) + adj_mat(*it, *kt) + adj_mat(*kt, *lt);
            auto ld = adj_mat(*lt, *jt) + adj_mat(*lt, *kt) + adj_mat(*it, *lt);
            if (id == 2 && jd == 2 && kd == 2 && ld == 2) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_five_almost_cliques(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;

  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            auto id = adj_mat(*it, *jt) + adj_mat(*it, *kt) + adj_mat(*it, *lt);
            auto jd = adj_mat(*it, *jt) + adj_mat(*jt, *kt) + adj_mat(*jt, *lt);
            auto kd = adj_mat(*kt, *jt) + adj_mat(*it, *kt) + adj_mat(*kt, *lt);
            auto ld = adj_mat(*lt, *jt) + adj_mat(*lt, *kt) + adj_mat(*it, *lt);
            if (id + jd + kd + ld == 10) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 3;
}

uint64_t count_five_cliques(const Graph &graph, AdjMatrixFn const &adj_mat) {
  uint64_t result = 0;

  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!adj_mat(*it, *jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!adj_mat(*it, *kt) || !adj_mat(*jt, *kt)) {
            continue;
          }

          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (!adj_mat(*it, *lt) || !adj_mat(*jt, *lt) || !adj_mat(*kt, *lt)) {
              continue;
            }

            result++;
          }
        }
      }
    }
  }
  return result / 3;
}

std::array<uint64_t, GRAPHLETS.size()> count_small_graphlets(const Graph &graph) {
  std::array<uint64_t, GRAPHLETS.size()> result;
  result.fill(0);

  std::vector<std::vector<bool>> adj_matrix(graph.size(), std::vector<bool>(graph.size(), false));
  for (size_t i = 0; i < graph.size(); i++) {
    for (auto &&v : graph[i]) {
      adj_matrix[i][v] = true;
    }
  }

  std::function edge_exists = [&](const unsigned &u, const unsigned &v) -> bool {
    return adj_matrix[u][v];
  };

  for (size_t i = 0; i < result.size(); i++) {
    result[i] = count_graphlets(graph, GRAPHLETS[i], &edge_exists);
  }

  return result;
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
  next_perm : { }
  } while (std::next_permutation(combination.begin(), combination.end()));
  return false;
}

uint64_t count_graphlets_naive(const Graph &graph, const Graph &graphlet) {
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

  uint64_t result = 0;
  do {
    if (is_isomorphic(graph, graphlet, combination)) {
      result++;
    }
  } while (next_combination());
  return result;
}
