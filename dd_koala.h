#pragma once

#include <graph/graph.h>

#include <set>
#include <vector>

inline int get_graph_size(const Koala::Graph<int, int> &G) {
  return G.getVertNo();
}

inline int get_index(
    const Koala::Graph<int, int>&, const Koala::Graph<int, int>::PVertex &v) {
  return v->getInfo();
}

inline int get_index(
    const Koala::Graph<int, int>::PVertex &u, const Koala::Graph<int, int>::PVertex &v,
    const int &n) {
  return u->getInfo() * n + v->getInfo();
}

inline void set_index(
    Koala::Graph<int, int>&, Koala::Graph<int, int>::PVertex &v, const int &value) {
  v->setInfo(value);
}

inline int get_degree(
    const Koala::Graph<int, int> &G, const Koala::Graph<int, int>::PVertex &v) {
  return G.deg(v);
}

inline void add_edge(
    Koala::Graph<int, int> &G, const Koala::Graph<int, int>::PVertex &v,
    const Koala::Graph<int, int>::PVertex &u) {
  G.addEdge(v, u);
}

inline bool check_edge(
    Koala::Graph<int, int> &G, const Koala::Graph<int, int>::PVertex &v,
    const Koala::Graph<int, int>::PVertex &u) {
  return G.getEdge(v, u) != NULL;
}

inline Koala::Graph<int, int>::PVertex add_vertex(Koala::Graph<int, int> &G, const int &value) {
  return G.addVert(value);
}

inline void move_vertex(
    Koala::Graph<int, int> &G, Koala::Graph<int, int> &H,
    Koala::Graph<int, int>::PVertex &v) {
  G.move(H, v);
}

inline void delete_vertex(
    Koala::Graph<int, int> &G, Koala::Graph<int, int>::PVertex &v) {
  G.delVert(v);
}

std::vector<Koala::Graph<int, int>::PVertex> get_vertices(const Koala::Graph<int, int> &G) {
  std::vector<Koala::Graph<int, int>::PVertex> V(get_graph_size(G));
  G.getVerts(V.begin());
  return V;
}

inline std::set<Koala::Graph<int, int>::PVertex> get_neighbors(
    const Koala::Graph<int, int> &G, const Koala::Graph<int, int>::PVertex &v) {
  std::set<Koala::Graph<int, int>::PVertex> S = G.getNeighSet(v);
  return S;
}
