#pragma once

#include "./dd_header.h"

typedef std::vector<std::set<unsigned>> Graph;

class DataObject {
 public:
  double open_triangles = 0, triangles = 0, average_degree = 0, average_degree_squared = 0;
  int no_vertices = 0, low_degree = 0, high_degree = 0;

  void print() const {
    printf("%5d vertices, %5d min degree, %5d max degree\n", no_vertices, low_degree, high_degree);
    printf("%10.0lf open triangles, %10.0lf triangles\n", open_triangles, triangles);
    printf(
        "%10.3lf average degree, %10.3lf average degree squared\n",
        average_degree, average_degree_squared);
  }

  std::string to_string() const {
    return "average degree: " + std::to_string(average_degree) + ", "
        + "open triangles: " + std::to_string(open_triangles) + ", "
        + "triangles: " + std::to_string(triangles);
  }
};

void degree_distribution(const Graph &G, DataObject &data) {
  data.no_vertices = G.size(), data.low_degree = G.size(), data.high_degree = 0;
  for (auto v : G) {
    data.low_degree = std::min(data.low_degree, static_cast<int>(v.size()));
    data.high_degree = std::max(data.high_degree, static_cast<int>(v.size()));
    data.average_degree += v.size();
    data.average_degree_squared += (v.size() * static_cast<double>(v.size()));
  }
  data.average_degree /= data.no_vertices;
  data.average_degree_squared /= data.no_vertices;
}

void count_triangles(const Graph &G, DataObject &data) {
  for (auto v : G) {
    data.open_triangles += v.size() * (v.size() - 1) / 2;

    for (auto it = v.begin(); it != v.end(); ++it) {
      for (auto jt = next(it, 1); jt != v.end(); ++jt) {
        if (G[*it].count(*jt)) {
          data.triangles++;
        }
      }
    }
  }
  data.triangles /= 3;
}

DataObject get_params_for_graph(const Graph &G, bool verbose = false) {
  DataObject data;
  degree_distribution(G, data);
  count_triangles(G, data);

  if (verbose) {
    data.print();
  }
  return data;
}

DataObject get_params_for_synthetic_graph(
    const Graph &G0, const int &n, const Parameters &params) {
  Graph G = G0;
  generate_graph_simple(G, n, params);
  return get_params_for_graph(G);
}
