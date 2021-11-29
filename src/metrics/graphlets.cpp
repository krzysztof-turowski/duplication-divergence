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

uint64_t count_graphlets(const Graph &graph, const Graph &graphlet) {
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
  } else if (graphlet == GRAPHLETS[5]) {
    return count_triangles_with_antenna(graph);
  } else if (graphlet == GRAPHLETS[6]) {
    return count_four_almost_cliques(graph);
  } else if (graphlet == GRAPHLETS[7]) {
    return count_four_cliques(graph);
  } else if (graphlet == GRAPHLETS[8]) {
    return count_five_paths(graph);
  } else if (graphlet == GRAPHLETS[9]) {
    return count_three_stars_with_antenna(graph);
  } else if (graphlet == GRAPHLETS[10]) {
    return count_four_stars(graph);
  } else if (graphlet == GRAPHLETS[11]) {
    return count_triangles_with_two_antennas(graph);
  } else if (graphlet == GRAPHLETS[12]) {
    return count_triangles_with_long_antenna(graph);
  } else if (graphlet == GRAPHLETS[13]) {
    return count_four_stars_with_edge(graph);
  } else if (graphlet == GRAPHLETS[14]) {
    return count_polygons(graph);
  } else if (graphlet == GRAPHLETS[15]) {
    return count_squares_with_antenna(graph);
  } else if (graphlet == GRAPHLETS[16]) {
    return count_almost_four_cliques_with_antenna(graph);
  } else if (graphlet == GRAPHLETS[17]) {
    return count_bowties(graph);
  } else if (graphlet == GRAPHLETS[18]) {
    return count_almost_four_cliques_with_antenna_alt(graph);
  } else if (graphlet == GRAPHLETS[19]) {
    return count_clique_two_three(graph);
  } else if (graphlet == GRAPHLETS[20]) {
    return count_houses(graph);
  } else if (graphlet == GRAPHLETS[21]) {
    return count_clique_three_one_one(graph);
  } else if (graphlet == GRAPHLETS[22]) {
    return count_four_cliques_with_antenna(graph);
  } else if (graphlet == GRAPHLETS[23]) {
    return count_three_triangles(graph);
  } else if (graphlet == GRAPHLETS[24]) {
    return count_squares_with_three_cross(graph);
  } else if (graphlet == GRAPHLETS[25]) {
    return count_four_cliques_with_flag(graph);
  } else if (graphlet == GRAPHLETS[26]) {
    return count_subdivided_crosses(graph);
  } else if (graphlet == GRAPHLETS[27]) {
    return count_five_almost_cliques(graph);
  }
  return count_graphlets_naive(graph, graphlet);
}

uint64_t count_open_triangles(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          result++;
        }
      }
    }
  }
  return result;
}

uint64_t count_triangles(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          result++;
        }
      }
    }
  }
  return result / 3;
}

uint64_t count_four_paths(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
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

uint64_t count_three_stars(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!graph[*it].count(*jt) && !graph[*jt].count(*kt) && !graph[*kt].count(*it)) {
            result++;
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_squares(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
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

uint64_t count_triangles_with_antenna(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (graph[*it].count(*jt) + graph[*jt].count(*kt) + graph[*kt].count(*it) == 1) {
            result++;
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_almost_cliques(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (graph[*it].count(*jt) + graph[*jt].count(*kt) + graph[*kt].count(*it) == 2) {
            result++;
          }
        }
      }
    }
  }
  return result / 2;
}

uint64_t count_four_cliques(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.upper_bound(i); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (graph[*it].count(*jt) + graph[*jt].count(*kt) + graph[*kt].count(*it) == 3) {
            result++;
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_five_paths(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          std::vector<unsigned> ks;
          for (auto &&k : graph[*it]) {
            if (!v.count(k) && !graph[*jt].count(k)) {
              ks.push_back(k);
            }
          }
          for (auto &&l : graph[*jt]) {
            if (!v.count(l) && !graph[*it].count(l)) {
              for (auto &&k : ks) {
                if (!graph[l].count(k)) {
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

uint64_t count_three_stars_with_antenna(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!graph[*jt].count(*kt) && !graph[*kt].count(*it)) {
            for (auto &&l : graph[*it]) {
              if (!v.count(l) && !graph[*jt].count(l) && !graph[*kt].count(l)) {
                result++;
              }
            }

            for (auto &&l : graph[*jt]) {
              if (!v.count(l) && !graph[*it].count(l) && !graph[*kt].count(l)) {
                result++;
              }
            }

            for (auto &&l : graph[*kt]) {
              if (!v.count(l) && !graph[*it].count(l) && !graph[*jt].count(l)) {
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

uint64_t count_four_stars(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (graph[*jt].count(*kt) || graph[*it].count(*kt)) {
            continue;
          }

          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (!graph[*kt].count(*lt) && !graph[*lt].count(*it) && !graph[*jt].count(*lt)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_triangles_with_two_antennas(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          continue;
        }

        for (auto &&k : graph[*it]) {
          if (graph[*jt].count(k) != 0 || graph[i].count(k) != 0) {
            continue;
          }

          for (auto &&l : graph[*jt]) {
            if (graph[l].count(k) == 0 && graph[l].count(*it) == 0 && graph[l].count(i) == 0) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_triangles_with_long_antenna(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          continue;
        }

        for (auto kt = graph[*it].begin(); kt != graph[*it].end(); ++kt) {
          if (graph[*kt].count(*jt) || graph[*kt].count(i)) {
            continue;
          }
          for (auto lt = std::next(kt); lt != graph[*it].end(); ++lt) {
            if (!graph[*lt].count(*jt) && !graph[*lt].count(i) && graph[*lt].count(*kt)) {
              result++;
            }
          }
        }

        for (auto kt = graph[*jt].begin(); kt != graph[*jt].end(); ++kt) {
          if (graph[*kt].count(*it) || graph[*kt].count(i)) {
            continue;
          }
          for (auto lt = std::next(kt); lt != graph[*jt].end(); ++lt) {
            if (!graph[*lt].count(*it) && !graph[*lt].count(i) && graph[*lt].count(*kt)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_stars_with_edge(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (graph[*kt].count(*it) + graph[*kt].count(*jt) + graph[*kt].count(*lt)
                    + graph[*lt].count(*it) + graph[*lt].count(*jt) + graph[*it].count(*jt)
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

uint64_t count_polygons(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto &&j : v) {
      for (auto &&k : graph[j]) {
        if (v.count(k)) {
          continue;
        }
        for (auto &&l : graph[k]) {
          if (v.count(l) || graph[j].count(l)) {
            continue;
          }
          for (auto &&m : graph[l]) {
            if (v.count(m) && !graph[j].count(m) && !graph[k].count(m)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 10;
}

uint64_t count_squares_with_antenna(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && !v.count(k) && graph[*jt].count(k)) {
              for (auto &&l : v) {
                if (l != *it && l != *jt && l != k && !graph[l].count(*it) && !graph[l].count(*jt)
                    && !graph[l].count(k)) {
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

uint64_t count_almost_four_cliques_with_antenna(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && v.count(k) && graph[*jt].count(k)) {
              for (auto &&l : v) {
                if (l != *it && l != *jt && l != k && !graph[l].count(*it) && !graph[l].count(*jt)
                    && !graph[l].count(k)) {
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

uint64_t count_bowties(const Graph &graph) {
  uint64_t result = 0;
  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (graph[*it].count(*jt) && graph[*kt].count(*lt) && !graph[*it].count(*kt)
                && !graph[*it].count(*lt) && !graph[*jt].count(*kt) && !graph[*jt].count(*lt)) {
              result++;
            } else if (graph[*it].count(*kt) && graph[*jt].count(*lt) && !graph[*it].count(*jt)
                       && !graph[*it].count(*lt) && !graph[*kt].count(*jt)
                       && !graph[*kt].count(*lt)) {
              result++;
            } else if (graph[*it].count(*lt) && graph[*jt].count(*kt) && !graph[*it].count(*jt)
                       && !graph[*it].count(*kt) && !graph[*lt].count(*jt)
                       && !graph[*lt].count(*kt)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_almost_four_cliques_with_antenna_alt(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          for (auto &&k : graph[*it]) {
            if (i != k && !v.count(k) && graph[*jt].count(k)) {
              for (auto &&l : v) {
                if (l != *it && l != *jt && l != k && !graph[l].count(*it) && !graph[l].count(*jt)
                    && !graph[l].count(k)) {
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

uint64_t count_clique_two_three(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          continue;
        }

        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (graph[*it].count(*kt) || graph[*jt].count(*kt)) {
            continue;
          }

          for (auto &&l : graph[*it]) {
            if (i < l && l != *jt && l != *kt && graph[*jt].count(l) && graph[*kt].count(l)
                && !v.count(l)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_houses(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          continue;
        }

        for (auto &&k : graph[*it]) {
          if (graph[*jt].count(k) || graph[i].count(k)) {
            continue;
          }

          for (auto &&l : graph[*jt]) {
            if (graph[l].count(k) && !graph[l].count(*it) && !graph[l].count(i)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_clique_three_one_one(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (graph[*it].count(*jt)) {
          continue;
        }

        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (graph[*it].count(*kt) || graph[*jt].count(*kt)) {
            continue;
          }

          for (auto &&l : graph[*it]) {
            if (i < l && l != *jt && l != *kt && graph[*jt].count(l) && graph[*kt].count(l)
                && v.count(l)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_four_cliques_with_antenna(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.upper_bound(i); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!graph[*jt].count(*kt) || !graph[*kt].count(*it)) {
            continue;
          }

          for (auto &&l : v) {
            if (l != *it && l != *jt && l != *kt && !graph[*it].count(l) && !graph[*jt].count(l)
                && !graph[*kt].count(l)) {
              result++;
            }
          }
          for (auto &&l : graph[*it]) {
            if (l != i && l != *jt && l != *kt && !graph[i].count(l) && !graph[*jt].count(l)
                && !graph[*kt].count(l)) {
              result++;
            }
          }
          for (auto &&l : graph[*jt]) {
            if (l != *it && l != i && l != *kt && !graph[*it].count(l) && !graph[i].count(l)
                && !graph[*kt].count(l)) {
              result++;
            }
          }
          for (auto &&l : graph[*kt]) {
            if (l != *it && l != *jt && l != i && !graph[*it].count(l) && !graph[*jt].count(l)
                && !graph[i].count(l)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

uint64_t count_three_triangles(const Graph &graph) {
  uint64_t result = 0;

  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            auto id = graph[*it].count(*jt) + graph[*it].count(*kt) + graph[*it].count(*lt);
            auto jd = graph[*it].count(*jt) + graph[*jt].count(*kt) + graph[*jt].count(*lt);
            auto kd = graph[*kt].count(*jt) + graph[*it].count(*kt) + graph[*kt].count(*lt);
            auto ld = graph[*lt].count(*jt) + graph[*lt].count(*kt) + graph[*it].count(*lt);
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

uint64_t count_squares_with_three_cross(const Graph &graph) {
  uint64_t result = 0;

  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (i < *it && graph[*it].count(*jt) && graph[*it].count(*kt) && !graph[*jt].count(*kt)) {
            for (auto &&l : graph[*jt]) {
              if (!v.count(l) && l != *it && l != *kt && l != i && !graph[*it].count(l)
                  && graph[*kt].count(l)) {
                result++;
              }
            }
          }

          if (i < *jt && graph[*jt].count(*it) && graph[*jt].count(*kt) && !graph[*it].count(*kt)) {
            for (auto &&l : graph[*kt]) {
              if (!v.count(l) && l != *jt && l != *it && l != i && !graph[*jt].count(l)
                  && graph[*it].count(l)) {
                result++;
              }
            }
          }

          if (i < *kt && graph[*kt].count(*jt) && graph[*kt].count(*it) && !graph[*jt].count(*it)) {
            for (auto &&l : graph[*it]) {
              if (!v.count(l) && l != *kt && l != *jt && l != i && !graph[*kt].count(l)
                  && graph[*jt].count(l)) {
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

uint64_t count_four_cliques_with_flag(const Graph &graph) {
  uint64_t result = 0;
  for (size_t i = 0; i < graph.size(); i++) {
    const auto &v = graph[i];
    for (auto it = v.upper_bound(i); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!graph[*jt].count(*kt) || !graph[*kt].count(*it)) {
            continue;
          }

          for (auto &&l : v) {
            if (l != *it && l != *jt && l != *kt
                && (graph[*it].count(l) + graph[*jt].count(l) + graph[*kt].count(l) == 1)) {
              result++;
            }
          }
          for (auto &&l : graph[*it]) {
            if (l != i && l != *jt && l != *kt
                && (graph[i].count(l) + graph[*jt].count(l) + graph[*kt].count(l) == 1)) {
              result++;
            }
          }
          for (auto &&l : graph[*jt]) {
            if (l != *it && l != i && l != *kt
                && (graph[*it].count(l) + graph[i].count(l) + graph[*kt].count(l) == 1)) {
              result++;
            }
          }
          for (auto &&l : graph[*kt]) {
            if (l != *it && l != *jt && l != i
                && (graph[*it].count(l) + graph[*jt].count(l) + graph[i].count(l) == 1)) {
              result++;
            }
          }
        }
      }
    }
  }
  return result / 2;
}

uint64_t count_subdivided_crosses(const Graph &graph) {
  uint64_t result = 0;

  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            auto id = graph[*it].count(*jt) + graph[*it].count(*kt) + graph[*it].count(*lt);
            auto jd = graph[*it].count(*jt) + graph[*jt].count(*kt) + graph[*jt].count(*lt);
            auto kd = graph[*kt].count(*jt) + graph[*it].count(*kt) + graph[*kt].count(*lt);
            auto ld = graph[*lt].count(*jt) + graph[*lt].count(*kt) + graph[*it].count(*lt);
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

uint64_t count_five_almost_cliques(const Graph &graph) {
  uint64_t result = 0;

  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            auto id = graph[*it].count(*jt) + graph[*it].count(*kt) + graph[*it].count(*lt);
            auto jd = graph[*it].count(*jt) + graph[*jt].count(*kt) + graph[*jt].count(*lt);
            auto kd = graph[*kt].count(*jt) + graph[*it].count(*kt) + graph[*kt].count(*lt);
            auto ld = graph[*lt].count(*jt) + graph[*lt].count(*kt) + graph[*it].count(*lt);
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

uint64_t count_five_cliques(const Graph &graph) {
  uint64_t result = 0;

  for (const auto &v : graph) {
    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = std::next(it); jt != v.end(); ++jt) {
        if (!graph[*it].count(*jt)) {
          continue;
        }
        for (auto kt = std::next(jt); kt != v.end(); ++kt) {
          if (!graph[*it].count(*kt) || !graph[*jt].count(*kt)) {
            continue;
          }

          for (auto lt = std::next(kt); lt != v.end(); ++lt) {
            if (!graph[*it].count(*lt) || !graph[*jt].count(*lt) || !graph[*kt].count(*lt)) {
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
