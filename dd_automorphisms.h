#pragma once

extern "C" {
  #include "lib/nauty/nauty.h"
  #include "lib/nauty/nausparse.h"
  #include "lib/nauty/traces.h"
}

#include <stack>

double log_automorphisms_dense(const std::vector<std::set<int>> &G) {
  DYNALLSTAT(graph, g, g_size);
  DYNALLSTAT(int, lab, lab_size);
  DYNALLSTAT(int, ptn, ptn_size);
  DYNALLSTAT(int, orbits, orbits_size);
  statsblk stats;
  static DEFAULTOPTIONS_GRAPH(options);

  int n = G.size(), word = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE, word, n, NAUTYVERSIONID);

  DYNALLOC2(graph, g, g_size, word, n, "malloc");
  DYNALLOC1(int, lab, lab_size, n, "malloc");
  DYNALLOC1(int, ptn, ptn_size, n, "malloc");
  DYNALLOC1(int, orbits, orbits_size, n, "malloc");
  EMPTYGRAPH(g, word, n);
  for (int v = 0; v < n; v++) {
    for (auto u : G[v]) {
      ADDONEEDGE(g, v, u, word);
    }
  }
  densenauty(g, lab, ptn, orbits, &options, &stats, word, n, NULL);
  return log(stats.grpsize1) + stats.grpsize2 * log(10);
}

double log_automorphisms_sparse(const std::vector<std::set<int>> &G) {
  sparsegraph g;
  DYNALLSTAT(int, lab, lab_size);
  DYNALLSTAT(int, ptn, ptn_size);
  DYNALLSTAT(int, orbits, orbits_size);
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
  DYNALLOC1(int, lab, lab_size, n, "malloc");
  DYNALLOC1(int, ptn, ptn_size, n, "malloc");
  DYNALLOC1(int, orbits, orbits_size, n, "malloc");

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
  sparsenauty(&g, lab, ptn, orbits, &options, &stats, NULL);
  return log(stats.grpsize1) + stats.grpsize2 * log(10);
}

double log_automorphisms_traces(const std::vector<std::set<int>> &G) {
  DYNALLSTAT(int, lab, lab_size);
  DYNALLSTAT(int, ptn, ptn_size);
  DYNALLSTAT(int, orbits, orbits_size);
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
  DYNALLOC1(int, lab, lab_size, n, "malloc");
  DYNALLOC1(int, ptn, ptn_size, n, "malloc");
  DYNALLOC1(int, orbits, orbits_size, n, "malloc");

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
  Traces(&g, lab, ptn, orbits, &options, &stats, NULL);
  return log(stats.grpsize1) + stats.grpsize2 * log(10);
}

std::vector<int> connectivity(const std::vector<std::set<int>> &G) {
  std::vector<int> C(G.size()), sizes;
  for (int i = 0, count = 0; i < static_cast<int>(G.size()); i++) {
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

int isolated_nodes(const std::vector<std::set<int>> &G) {
  return std::count_if(G.begin(), G.end(), [](const std::set<int> &v){ return v.size() == 0; });
}
