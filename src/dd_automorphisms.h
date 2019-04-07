#pragma once

extern "C" {
  #include <nauty.h>
  #include <nausparse.h>
  #include <traces.h>
}

#include <set>
#include <stack>
#include <vector>

typedef std::vector<std::set<unsigned>> Graph;

double log_automorphisms_dense(const Graph &G) {
  statsblk stats;
  static DEFAULTOPTIONS_GRAPH(options);

  int n = G.size(), word = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE, word, n, NAUTYVERSIONID);

  std::vector<graph> g(word * n);
  std::vector<int> lab(n), ptn(n), orbits(n);
  EMPTYGRAPH(&g[0], word, n);
  for (int v = 0; v < n; v++) {
    for (auto u : G[v]) {
      ADDONEEDGE(&g[0], v, u, word);
    }
  }
  densenauty(&g[0], &lab[0], &ptn[0], &orbits[0], &options, &stats, word, n, NULL);
  return log(stats.grpsize1) + stats.grpsize2 * log(10);
}

double log_automorphisms_sparse(const Graph &G) {
  sparsegraph g;
  statsblk stats;
  static DEFAULTOPTIONS_SPARSEGRAPH(options);

  int n = G.size(), word = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE, word, n, NAUTYVERSIONID);

  int m = 0;
  for (int v = 0; v < n; v++) {
    m += G[v].size();
  }

  SG_INIT(g);
  SG_ALLOC(g, n, m, "malloc");
  std::vector<int> lab(n), ptn(n), orbits(n);

  g.nv = n;
  g.nde = m;
  int position = 0, index = 0;
  for (int v = 0; v < n; v++) {
    g.v[v] = position;
    g.d[v] = G[v].size();
    position += G[v].size();
    for (auto u : G[v]) {
      g.e[index] = u;
      index++;
    }
  }
  sparsenauty(&g, &lab[0], &ptn[0], &orbits[0], &options, &stats, NULL);
  return log(stats.grpsize1) + stats.grpsize2 * log(10);
}

double log_automorphisms_traces(const Graph &G) {
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;

  int n = G.size(), word = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE, word, n, NAUTYVERSIONID);

  int m = 0;
  for (int v = 0; v < n; v++) {
    m += G[v].size();
  }

  SG_DECL(g);
  SG_ALLOC(g, n, m, "malloc");
  std::vector<int> lab(n), ptn(n), orbits(n);

  g.nv = n;
  g.nde = m;
  int position = 0, index = 0;
  for (int v = 0; v < n; v++) {
    g.v[v] = position;
    g.d[v] = G[v].size();
    position += G[v].size();
    for (auto u : G[v]) {
      g.e[index] = u;
      index++;
    }
  }
  Traces(&g, &lab[0], &ptn[0], &orbits[0], &options, &stats, NULL);
  return log(stats.grpsize1) + stats.grpsize2 * log(10);
}

std::vector<int> connectivity(const Graph &G) {
  std::vector<int> C(G.size()), sizes;
  for (size_t i = 0, count = 0; i < G.size(); i++) {
    if (C[i] == 0) {
      int vertices = 0;
      std::stack<int> S;
      S.push(i), count++;
      while (!S.empty()) {
        int v = S.top();
        S.pop();
        if (C[v] != 0) {
          continue;
        }
        C[v] = count, vertices++;
        for (auto u : G[v]) {
          S.push(u);
        }
      }
      sizes.push_back(vertices);
    }
  }
  return sizes;
}

int isolated_nodes(const Graph &G) {
  return std::count_if(
      G.begin(), G.end(), [](const std::set<unsigned> &v){ return v.size() == 0; });
}
