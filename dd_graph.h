#pragma once

#include "./dd_header.h"

#if defined(koala)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-align"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wextra"
  #pragma GCC diagnostic ignored "-Wold-style-cast"
  #pragma GCC diagnostic ignored "-Wpedantic"
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wswitch-default"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wunused-variable"
  #if __GNUC__ >= 7
    #pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
  #endif
  #if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wmisleading-indentation"
  #endif
  #ifndef __clang__
    #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  #endif
  #include "./dd_koala.h"
  #pragma GCC diagnostic pop

  typedef Koala::Graph<int, int> Graph;
  typedef Koala::Graph<int, int>::PVertex Vertex;
#elif defined(snap)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wold-style-cast"
  #pragma GCC diagnostic ignored "-Wpedantic"
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wmisleading-indentation"
  #endif
  #include "./dd_snap.h"
  #pragma GCC diagnostic pop

  typedef TNodeNet<TInt> Graph;
  typedef int Vertex;
#elif defined(networkit)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wold-style-cast"
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #include "./dd_networkit.h"
  #pragma GCC diagnostic pop

  typedef NetworKit::Graph Graph;
  typedef NetworKit::node Vertex;
#endif

#include <cassert>
#include <limits>
#include <string>
#include <set>
#include <vector>

const double EPS = 10e-9;

Graph generate_seed(const int &n0, const double &p0) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  Graph G;
  std::vector<Vertex> V(n0);
  for (int i = 0; i < n0; i++) {
    V[i] = add_vertex(G, i);
  }
  for (int i = 0; i < n0; i++) {
    for (int j = i + 1; j < n0; j++) {
      if (edge_distribution(generator) <= p0) {
        add_edge(G, V[i], V[j]);
      }
    }
  }
  return G;
}

Graph generate_graph(Graph &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  std::vector<Vertex> V(get_vertices(G));
  V.resize(n);
  for (int i = get_graph_size(G); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    V[i] = add_vertex(G, i);
    std::set<Vertex> neighbors(get_neighbors(G, V[parent]));
    if (params.mode == Mode::PURE_DUPLICATION) {
      for (const auto &v : neighbors) {
        if (edge_distribution(generator) <= params.p) {
          add_edge(G, V[i], v);
        }
      }
    } else if (params.mode == Mode::PURE_DUPLICATION_CONNECTED) {
      while (true) {
        for (const auto &v : neighbors) {
          if (edge_distribution(generator) <= params.p) {
            add_edge(G, V[i], v);
          }
        }
        if (get_degree(G, V[i]) > 0) {
          break;
        }
        neighbors = get_neighbors(G, V[parent_distribution(generator)]);
      }
    } else if (params.mode == Mode::CHUNG_LU) {
      for (const auto &v : neighbors) {
        if (edge_distribution(generator) <= params.p) {
          add_edge(G, V[i], v);
        }
      }
      if (edge_distribution(generator) <= params.q) {
        add_edge(G, V[i], V[parent]);
      }
    } else if (params.mode == Mode::PASTOR_SATORRAS) {
      for (int j = 0; j < i; j++) {
        if (neighbors.count(V[j])) {
          if (edge_distribution(generator) <= params.p) {
            add_edge(G, V[i], V[j]);
          }
        } else {
          if (edge_distribution(generator) <= params.r / i) {
            add_edge(G, V[i], V[j]);
          }
        }
      }
    } else {
      throw std::invalid_argument("Invalid mode: " + params.to_string());
    }
  }
  return G;
}

Graph read_graph(const std::string &graph_name) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  Graph G;
  std::vector<Vertex> V;
  int u, v;
  while (graph_file >> u >> v) {
    if (v >= get_graph_size(G)) {
      for (int i = get_graph_size(G); i <= v; i++) {
        V.push_back(add_vertex(G, i));
      }
    }
    if (u != v) {
      add_edge(G, V[u], V[v]);
    }
  }
  graph_file.close();
  return G;
}

class NeighborhoodStructure {
 protected:
  template <typename T>
  struct counting_iterator {
    size_t count;
    T dummy;

    counting_iterator() : count(0) { }
    counting_iterator& operator++() { ++count; return *this; }
    counting_iterator operator++(int) { ++count; return *this; }
    T& operator*() { return dummy; }
  };

 public:
  virtual int common_neighbors(const Vertex &v, const Vertex &u) const = 0;
  virtual void remove_vertex(const std::set<Vertex> &neighbors) = 0;
  virtual void restore_vertex(const std::set<Vertex> &neighbors) = 0;
  virtual bool verify(const Graph &G) const = 0;
};

class NoNeighborhoodStructure : public NeighborhoodStructure {
 private:
  const Graph &G;

 public:
  explicit NoNeighborhoodStructure(const Graph &H) : G(H) { }

  int common_neighbors(const Vertex &v, const Vertex &u) const {
    auto N_v(get_neighbors(G, v));
    auto N_u(get_neighbors(G, u));
    return set_intersection(
        N_v.begin(), N_v.end(), N_u.begin(), N_u.end(), counting_iterator<Vertex>()).count;
  }

  void remove_vertex(const std::set<Vertex>&) { }

  void restore_vertex(const std::set<Vertex>&) { }

  bool verify(const Graph &H) const { return &G == &H; }
};

class CompleteNeighborhoodStructure : public NeighborhoodStructure {
 private:
  int n;
  std::vector<int> V;

 public:
  explicit CompleteNeighborhoodStructure(
      const Graph &G) : n(get_graph_size(G)), V(get_graph_size(G) * get_graph_size(G)) {
    auto vertices = get_vertices(G);
    #pragma omp parallel for
    for (std::size_t v_i = 0; v_i < vertices.size(); v_i++) {
      auto N_v(get_neighbors(G, vertices[v_i]));
      for (std::size_t u_i = 0; u_i < vertices.size(); u_i++) {
        auto N_u(get_neighbors(G, vertices[u_i]));
        V[get_index(vertices[v_i], vertices[u_i], n)] =
            set_intersection(
                N_v.begin(), N_v.end(), N_u.begin(), N_u.end(),
                counting_iterator<Vertex>()).count;
      }
    }
  }

  explicit CompleteNeighborhoodStructure(
      const CompleteNeighborhoodStructure &other) : n(other.n), V(other.V) { }

  int common_neighbors(const Vertex &v, const Vertex &u) const {
    return V[get_index(v, u, n)];
  }

  void remove_vertex(const std::set<Vertex> &neighbors) {
    for (const auto &w : neighbors) {
      for (const auto &u : neighbors) {
        --V[get_index(w, u, n)];
      }
    }
  }

  void restore_vertex(const std::set<Vertex> &neighbors) {
    for (const auto &w : neighbors) {
      for (const auto &u : neighbors) {
        ++V[get_index(w, u, n)];
      }
    }
  }

  bool verify(const Graph &G) const {
    auto vertices = get_vertices(G);
    for (const auto &v : vertices) {
      std::set<Vertex> N_v(get_neighbors(G, v));
      for (const auto &u : vertices) {
        std::set<Vertex> N_u(get_neighbors(G, u));
        int common_from_graph =
            set_intersection(
                N_v.begin(), N_v.end(), N_u.begin(), N_u.end(),
                counting_iterator<Vertex>()).count;
        int common_from_struct = V[get_index(v, u, n)];
        if (common_from_graph != common_from_struct) {
          return false;
        }
      }
    }
    return true;
  }
};

inline long double add_exp_log(const long double &x, const long double &y) {
  if (x == -std::numeric_limits<long double>::infinity()) {
    return y;
  }
  if (y == -std::numeric_limits<long double>::infinity()) {
    return x;
  }
  return x > y ? x + log2l(1.0L + exp2l(y - x)) : y + log2l(1.0L + exp2l(x - y));
}

bool is_feasible(
    const Graph &G, const Parameters &params, const Vertex &v, const Vertex &u,
    const NeighborhoodStructure &aux) {
  switch (params.mode) {
    case Mode::PURE_DUPLICATION: {
      bool uv = check_edge(G, u, v);
      int both = aux.common_neighbors(v, u), only_v = get_degree(G, v) - both - uv;
      if (!uv && only_v == 0) {
        return true;
      }
      return false;
    }
    case Mode::PASTOR_SATORRAS:
      return true;
    default:
      throw std::invalid_argument("Invalid mode: " + params.to_string());
  }
}

bool is_feasible(
    const Graph &G, const Parameters &params, const Vertex &v,
    const NeighborhoodStructure &aux) {
  switch (params.mode) {
    case Mode::PURE_DUPLICATION: {
      std::vector<Vertex> V(get_vertices(G));
      for (const auto &u : V) {
        if (u != v) {
          if (is_feasible(G, params, v, u, aux)) {
            return true;
          }
        }
      }
      return false;
    }
    case Mode::PASTOR_SATORRAS:
      return true;
    default:
      throw std::invalid_argument("Invalid mode: " + params.to_string());
  }
}

long double get_log_transition_probability(
    const Graph &G, const Parameters &params, const Vertex &v, const Vertex &u,
    const NeighborhoodStructure &aux) {
  if (!is_feasible(G, params, v, u, aux)) {
    return -std::numeric_limits<long double>::infinity();
  }
  bool uv = check_edge(G, u, v);
  int both = aux.common_neighbors(v, u), only_v = get_degree(G, v) - both,
      only_u = get_degree(G, u) - both - uv,
      none = (get_graph_size(G) - 1) + both - only_u - only_v;
  long double p(params.p), r(params.r);
  switch (params.mode) {
    case Mode::PURE_DUPLICATION:
      assert(!(fabsl(p) < EPS && both > 0) && !(fabsl(1 - p) < EPS && only_u > 0));
      return both * log2l(p) + only_u * log2l(1 - p) - log2l(get_graph_size(G) - 1);
    case Mode::PASTOR_SATORRAS:
      assert(!(fabsl(p) < EPS && both > 0) && !(fabsl(1 - p) < EPS && only_u > 0)
          && !(fabsl(r) < EPS && only_v > 0)
          && !(fabsl((get_graph_size(G) - 1) - r) < EPS && none > 0));
      return both * log2l(p) + only_u * log2l(1 - p) + only_v * log2l(r)
          + none * log2l(get_graph_size(G) - 1 - r)
          - (only_v + none + 1) * log2l(get_graph_size(G) - 1);
    default:
      throw std::invalid_argument("Invalid mode: " + params.to_string());
  }
}

long double get_log_transition_probability(
    const Graph &G, const Parameters &params,
    const Vertex &v, const NeighborhoodStructure &aux) {
  long double p_v = -std::numeric_limits<long double>::infinity();
  std::vector<Vertex> V(get_vertices(G));
  for (const auto &u : V) {
    if (u != v) {
      p_v = add_exp_log(p_v, get_log_transition_probability(G, params, v, u, aux));
    }
  }
  return p_v;
}

std::vector<long double> get_transition_probability(
    const Graph &G, const Parameters &params, const NeighborhoodStructure &aux) {
  std::vector<long double> out;
  std::vector<Vertex> V(get_vertices(G));
  for (const auto &v : V) {
    out.push_back(exp2l(get_log_transition_probability(G, params, v, aux)));
  }
  return out;
}

long double get_discard_score(
    const Graph &G, const Parameters &params, const Vertex &v, const Vertex &u,
    const NeighborhoodStructure &aux) {
  bool uv = check_edge(G, u, v);
  int both = aux.common_neighbors(v, u), only_v = get_degree(G, v) - both - uv,
      only_u = get_degree(G, u) - both - uv;
  switch (params.mode) {
    case Mode::PURE_DUPLICATION:
      return !uv && only_v == 0 ? expl(-4 * only_u) : 0.0L;
    case Mode::PASTOR_SATORRAS:
      return expl(-4 * only_u - 8 * only_v);
    default:
      throw std::invalid_argument("Invalid mode: " + params.to_string());
  }
}

long double get_discard_score(
    const Graph &G, const Parameters &params, const Vertex &v,
    const NeighborhoodStructure &aux) {
  long double out = 0.0L;
  std::vector<Vertex> V(get_vertices(G));
  for (const auto &u : V) {
    out += get_discard_score(G, params, v, u, aux);
  }
  return out;
}
