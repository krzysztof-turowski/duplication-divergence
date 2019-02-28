#pragma once

#include <set>
#include <vector>

#include "./lib/snap/snap-core/Snap.h"

inline int get_graph_size(const TNodeNet<TInt> &G) {
  return G.GetNodes();
}

inline int get_index(const int &v) {
  return v;
}

inline int get_index(const int &u, const int &v, const int &n) {
  return v * n + u;
}

inline void set_index(TNodeNet<TInt> &G, int &v, const int &value) {
  G.GetNI(v).GetDat() = value;
}

inline int get_degree(TNodeNet<TInt> &G, const int &v) {
  return G.GetNI(v).GetDeg();
}

inline void add_edge(TNodeNet<TInt> &G, const int &v, const int &u) {
  G.AddEdge(v, u);
}

inline bool check_edge(TNodeNet<TInt> &G, const int &v, const int &u) {
  return G.IsEdge(v, u);
}

inline int add_vertex(TNodeNet<TInt> &G, const int &value) {
  auto v = G.AddNode(value);
  set_index(G, v, value);
  return v;
}

inline void move_vertex(TNodeNet<TInt> &G, TNodeNet<TInt> &H, int &v) {
  H.DelNode(v), G.AddNode(v);
}

inline void delete_vertex(TNodeNet<TInt> &G, int &v) {
  G.DelNode(v);
}

std::vector<int> get_vertices(const TNodeNet<TInt> &G) {
  std::vector<int> V;
  for (auto it = G.BegNI(); it < G.EndNI(); it++) {
    V.push_back(it.GetId());
  }
  return V;
}

inline std::set<int> get_neighbors(const TNodeNet<TInt> &G, const int &v) {
  std::set<int> S;
  auto vertex = G.GetNI(v);
  for (int e = 0; e < vertex.GetDeg(); e++) {
    S.insert(vertex.GetNbrNId(e));
  }
  return S;
}
