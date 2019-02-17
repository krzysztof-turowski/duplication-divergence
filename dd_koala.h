#pragma once

#include "dd_header.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-default"
#include "lib/koala/graph/graph.h"
#pragma GCC diagnostic pop

Koala::Graph<int, int> generate_seed_koala(const int &n0, const double &p0) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  Koala::Graph<int, int> G;
  std::vector<Koala::Graph<int, int>::PVertex> V;
  for (int i = 0; i < n0; i++) {
    V.push_back(G.addVert(i));
  }
  for (int i = 0; i < n0; i++) {
    for (int j = i + 1; j < n0; j++) {
      if (edge_distribution(generator) <= p0) {
        G.addEdge(V[i], V[j]);
      }
    }
  }
  return G;
}

Koala::Graph<int, int> generate_graph_koala(Koala::Graph<int, int> &G, const int &n, const Parameters &params) {
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> edge_distribution(0.0, 1.0);

  std::vector<Koala::Graph<int, int>::PVertex> V(n);
  G.getVerts(V.begin());
  for (int i = G.getVertNo(); i < n; i++) {
    std::uniform_int_distribution<int> parent_distribution(0, i - 1);
    int parent = parent_distribution(generator);
    V[i] = G.addVert(i);
    std::set<Koala::Graph<int, int>::PVertex> neighbors = G.getNeighSet(V[parent]);
    if (params.mode == Mode::PURE_DUPLICATION) {
      for (auto v : neighbors) {
        if (edge_distribution(generator) <= params.p) {
          G.addEdge(V[i], v);
        }
      }
    }
    else if (params.mode == Mode::PURE_DUPLICATION_CONNECTED) {
      while(true) {
        for (auto v : neighbors) {
          if (edge_distribution(generator) <= params.p) {
            G.addEdge(V[i], v);
          }
        }
        if (G.deg(V[i]) > 0) {
          break;
        }
        neighbors = G.getNeighSet(V[parent_distribution(generator)]);
      }
    }
    else if (params.mode == Mode::CHUNG_LU) {
      for (auto v : neighbors) {
        if (edge_distribution(generator) <= params.p) {
          G.addEdge(V[i], v);
        }
      }
      if (edge_distribution(generator) <= params.q) {
        G.addEdge(V[i], V[parent]);
      }
    }
    else if (params.mode == Mode::PASTOR_SATORRAS) {
      for (int j = 0; j < i; j++) {
        if (neighbors.count(V[j])) {
          if (edge_distribution(generator) <= params.p) {
            G.addEdge(V[i], V[j]);
          }
        }
        else {
          if (edge_distribution(generator) <= params.r / i) {
            G.addEdge(V[i], V[j]);
          }
        }
      }
    }
    else {
      assert(0);
    }
  }
  return G;
}

Koala::Graph<int, int> read_graph_koala(const std::string &graph_name) {
  std::ifstream graph_file(graph_name);
  if (graph_file.fail()) {
    throw std::invalid_argument("Missing " + graph_name + " file");
  }
  Koala::Graph<int, int> G;
  std::vector<Koala::Graph<int, int>::PVertex> V;
  int u, v;
  while (!graph_file.eof())
  {
    graph_file >> u >> v;
    if (v >= G.getVertNo()) {
      for (int i = G.getVertNo(); i <= v; i++) {
        V.push_back(G.addVert(i));
      }
    }
    if (u != v) {
      G.addEdge(V[u], V[v]);
    }
  }
  graph_file.close();
  return G;
}

template <typename T>
struct counting_iterator
{
    size_t count;
    T dummy;

    counting_iterator() : count(0) { }
    counting_iterator& operator++() { ++count; return *this; }
    T& operator*() { return dummy; }
};

class NeighborhoodStructure {
private:
  int n;
  std::vector<int> V;

  inline int get_index(const Koala::Graph<int, int>::PVertex &v, const Koala::Graph<int, int>::PVertex &u) const {
    return u->getInfo() * n + v->getInfo();
  }

public:
  NeighborhoodStructure(const Koala::Graph<int, int> &G) : n(G.getVertNo()) {
    V.resize(n * n);
    for (auto v = G.getVert(); v; v = G.getVertNext(v)) {
      std::set<Koala::Graph<int, int>::PVertex> N_v = G.getNeighSet(v);
      for (auto u = G.getVert(); u; u = G.getVertNext(u)) {
        std::set<Koala::Graph<int, int>::PVertex> N_u = G.getNeighSet(u);
        V[get_index(v, u)] = set_intersection(N_v.begin(), N_v.end(), N_u.begin(), N_u.end(), counting_iterator<Koala::Graph<int, int>::PVertex>()).count;
      }
    }
  }

  int common_neighbors(const Koala::Graph<int, int>::PVertex &v, const Koala::Graph<int, int>::PVertex &u) const {
    return V[get_index(v, u)];
  }

  void remove_vertex(const Koala::Graph<int, int>::PVertex &v, const std::set<Koala::Graph<int, int>::PVertex> neighbors) {
    for (auto w : neighbors) {
      for (auto u : neighbors) {
        --V[get_index(w, u)];
      }
    }
  }

  void restore_vertex(const Koala::Graph<int, int>::PVertex &v, const std::set<Koala::Graph<int, int>::PVertex> neighbors) {
    for (auto w : neighbors) {
      for (auto u : neighbors) {
        ++V[get_index(w, u)];
      }
    }
  }
};

double get_transition_probability(
    const Koala::Graph<int, int> &G, const Parameters &params, const Koala::Graph<int, int>::PVertex &v, const Koala::Graph<int, int>::PVertex &u,
    const NeighborhoodStructure &aux) {
  switch (params.mode) {
    case Mode::PURE_DUPLICATION:
      if (G.deg(v) <= G.deg(u)) {
        int both = aux.common_neighbors(v, u);
        return G.deg(v) - both == 0 ? pow(params.p, both) * pow(1 - params.p, G.deg(u) - both) / (G.getVertNo() - 1) : 0.0;
      }
      return 0.0;
    case Mode::PASTOR_SATORRAS: {
        int both = aux.common_neighbors(v, u), only_v = G.deg(v) - both, only_u = G.deg(u) - both, none = G.getVertNo() + both - only_u - only_v;
        return pow(params.p, both) * pow(params.r / (G.getVertNo() - 1), only_v) * pow(1 - params.p, only_u) * pow(1 - (params.r / (G.getVertNo() - 1)), none)
            / (G.getVertNo() - 1);
      }
    default:
      throw std::invalid_argument("Invalid mode: " + params.to_string());
  }
}

double get_transition_probability(
    const Koala::Graph<int, int> &G, const Parameters &params, const Koala::Graph<int, int>::PVertex &v,
    const NeighborhoodStructure &aux) {
  double p_v = 0;
  for (auto u = G.getVert(); u; u = G.getVertNext(u)) {
    if (u != v) {
      p_v += get_transition_probability(G, params, v, u, aux);
    }
  }
  return p_v;
}

std::vector<double> get_transition_probability(const Koala::Graph<int, int> &G, const Parameters &params, const NeighborhoodStructure &aux) {
  std::vector<double> out;
  for (Koala::Graph<int, int>::PVertex v = G.getVert(); v; v = G.getVertNext(v)) {
    out.push_back(get_transition_probability(G, params, v, aux));
  }
  return out;
}
