#pragma once

#include <graph/Graph.h>

#include <set>
#include <vector>

// TODO(kturowski): implement labels as a pair (NetworKit::Graph, map<NetworKit::node, int>)

inline int get_graph_size(const NetworKit::Graph &G) {
  return G.numberOfNodes();
}

inline int get_index(const NetworKit::Graph &G, const NetworKit::node &v) {
  throw std::logic_error("Not yet implemented");
}

inline int get_index(const NetworKit::node &u, const NetworKit::node &v, const int &n) {
  throw v * n + u;
}

inline void set_index(NetworKit::Graph &G, NetworKit::node &v, const int &value) {
  throw std::logic_error("Not yet implemented");
}

inline int get_degree(const NetworKit::Graph &G, const NetworKit::node &v) {
  return G.degree(v);
}

inline void add_edge(NetworKit::Graph &G, const NetworKit::node &v, const NetworKit::node &u) {
  G.addEdge(v, u);
}

inline bool check_edge(
    const NetworKit::Graph &G, const NetworKit::node &v, const NetworKit::node &u) {
  return G.hasEdge(v, u);
}

inline NetworKit::node add_vertex(NetworKit::Graph &G, const int &value) {
  return G.addNode();
}

inline void move_vertex(NetworKit::Graph &G, NetworKit::Graph &H, NetworKit::node &v) {
  if (G.isEmpty()) {
    H.removeNode(v), G.addNode();
  } else {
    G.restoreNode(v), H.removeNode(H.randomNode());
  }
}

inline void delete_vertex(NetworKit::Graph &G, NetworKit::node &v) {
  G.removeNode(v);
}

std::vector<NetworKit::node> get_vertices(const NetworKit::Graph &G) {
  return G.nodes();
}

inline std::set<NetworKit::node> get_neighbors(
    const NetworKit::Graph &G, const NetworKit::node &v) {
  std::vector<NetworKit::node> V(G.neighbors(v));
  std::set<NetworKit::node> S(V.begin(), V.end());
  return S;
}
