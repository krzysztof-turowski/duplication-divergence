#pragma once

#include <set>
#include <vector>

#include "./lib/snap/snap-core/Snap.h"

inline int get_graph_size(const TNodeNet<TInt> &G) {
  return G.GetNodes();
}

inline int get_index(const TNodeNet<TInt>::TNodeI &v) {
  return v.GetDat();
}

inline int get_index(
    const TNodeNet<TInt>::TNodeI &u, const TNodeNet<TInt>::TNodeI &v, const int &n) {
  return v.GetDat() * n + u.GetDat();
}

inline void set_index(TNodeNet<TInt>::TNodeI &v, const int &value) {
  v.GetDat() = value;
}

inline int get_degree(TNodeNet<TInt> &G, const TNodeNet<TInt>::TNodeI &v) {
  return v.GetDeg();
}

inline void add_edge(
    TNodeNet<TInt> &G, const TNodeNet<TInt>::TNodeI &v, const TNodeNet<TInt>::TNodeI &u) {
  G.AddEdge(v.GetId(), u.GetId());
}

inline bool check_edge(
    TNodeNet<TInt> &G, const TNodeNet<TInt>::TNodeI &v, const TNodeNet<TInt>::TNodeI &u) {
  return G.IsEdge(v.GetId(), u.GetId());
}

inline TNodeNet<TInt>::TNodeI add_vertex(TNodeNet<TInt> &G, const int &value) {
  auto v = G.GetNI(G.AddNode());
  set_index(v, value);
  return v;
}

inline void move_vertex(TNodeNet<TInt> &G, TNodeNet<TInt> &H, TNodeNet<TInt>::TNodeI &v) {
  if (G.Empty()) {
    H.DelNode(v.GetId());
  } else {
    throw std::logic_error("Restore not yet implemented");
  }
}

inline void delete_vertex(TNodeNet<TInt> &G, TNodeNet<TInt>::TNodeI &v) {
  G.DelNode(v.GetId());
}

std::vector<TNodeNet<TInt>::TNodeI> get_vertices(const TNodeNet<TInt> &G) {
  std::vector<TNodeNet<TInt>::TNodeI> V;
  for (auto it = G.BegNI(); it < G.EndNI(); it++) {
    V.push_back(it);
  }
  return V;
}

inline std::set<TNodeNet<TInt>::TNodeI> get_neighbors(
    const TNodeNet<TInt> &G, const TNodeNet<TInt>::TNodeI &v) {
  std::set<TNodeNet<TInt>::TNodeI> S;
  for (int e = 0; e < v.GetOutDeg(); e++) {
    S.insert(G.GetNI(v.GetOutNId(e)));
  }
  return S;
}
